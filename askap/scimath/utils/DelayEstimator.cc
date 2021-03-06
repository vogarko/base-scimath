/// @file
/// 
/// @brief estimate delay from the complex spectrum 
/// @details This class implements a simple algorithm to estimate delay from a complex spectrum
/// by unwrapping the phase and fitting a straight line into the phase slope. This code was
/// originally located in the data monitor of the software correlator, but was moved here so
/// we can reuse it in the CP ingest pipeline.  
///
///
/// @copyright (c) 2007 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#include "askap/scimath/utils/DelayEstimator.h"
#include "askap/scimath/utils/PhaseUnwrapper.h"
#include "askap/scimath/fft/FFTWrapper.h"
#include <casacore/casa/BasicSL/Constants.h>
#include <askap/askap/AskapError.h>

#include <vector>
#include <cmath>

namespace askap {

namespace scimath {

/// @brief construct estimator for a given spectral resolution
/// @param[in] resolution the spectral resolution in Hz
DelayEstimator::DelayEstimator(const double resolution) : itsResolution(resolution), itsQuality(0.) {}
   
/// @brief estimate delay for a given spectrum
/// @param[in] vis (visibility) spectrum
/// @return delay in seconds
double DelayEstimator::getDelay(const casacore::Vector<casacore::Complex> &vis) const
{
   ASKAPASSERT(itsResolution != 0.);
   ASKAPASSERT(vis.nelements() > 1);
   std::vector<float> phases(vis.nelements());
   const float threshold = 3 * casacore::C::pi / 2;

   // unambiguate phases
   PhaseUnwrapper<float> phu(threshold);
   for (size_t chan=0; chan<phases.size(); ++chan) {
        const float curPhase = arg(vis[chan]);
        if (std::isnan(curPhase)) {
            phases[chan] = curPhase;
        } else {
            phases[chan] = phu(arg(vis[chan]));
        }
   }
   
   // do LSF into phase vs. channel
   double sx = 0., sy = 0., sx2 = 0., sxy = 0., sy2 = 0.;
   // could've combined two loops, but keep it easy for now
   for (size_t chan=0; chan < phases.size(); ++chan) {
        if (std::isnan(phases[chan])) {
            continue;
        }
        const double phase = double(phases[chan]); 
        sx += double(chan);
        sx2 += double(chan) * double(chan);
        sy += phase;
        sxy += double(chan) * phase;
        sy2 += phase * phase;
   }
   sx /= double(phases.size());
   sy /= double(phases.size());
   sx2 /= double(phases.size());
   sxy /= double(phases.size());
   sy2 /= double(phases.size());
   const double coeffNumerator = sxy - sx * sy;
   const double chanVariance = sx2 - sx * sx;
   ASKAPDEBUGASSERT(chanVariance != 0.);
   const double coeff =  coeffNumerator / chanVariance;
   // for this method quality is determined by the chi-squared. It is found better than 
   // the correlation coefficient, because the latter could be quite close to zero if the delay is small but
   // bandpass shape cannot be approximated by a line. the absolute value of the correlation coefficient
   const double phaseVariance = sy2 - sy * sy;
   if (phaseVariance > 0.) {
       const double chi2 = phaseVariance - 2 * coeff * coeffNumerator + casacore::square(coeff) * chanVariance;
       // abs value of correlation coefficient: fabs(coeffNumerator / sqrt(chanVariance * phaseVariance));
       // atan2 is just a convenient non-linear function to map an inverse into [0,1] range and avoid issues
       // with infinities
       itsQuality = atan2(1., chi2)  * casacore::C::_2_pi;
       // calculate delay based on the fitted slope
       return coeff / 2. / casacore::C::pi / itsResolution;
   }
   // degenerate case - zero delay
   itsQuality = 1.;
   return 0.;
}

/// @brief estimate delay for a given spectrum via FFT
/// @details This method works well in the case of multiple harmonics
/// present. However, it only gives a rough estimate
/// @param[in] vis (visibility) spectrum
/// @return delay in seconds
double DelayEstimator::getDelayWithFFT(const casacore::Vector<casacore::Complex> &vis) const
{
  ASKAPASSERT(itsResolution != 0.);
   
  // create a copy explicitly due to reference semantics of casa arrays
  casacore::Vector<casacore::Complex> lags(vis.copy());
  scimath::fft(lags, true);
  // search for a peak lag
  casacore::uInt peakLagChan = lags.nelements(); 
  float peakAmp = -1.;
  float meanAmp = 0.;
  for (casacore::uInt chan = 0; chan < lags.nelements(); ++chan) {
       const float curAmp = abs(lags[chan]);
       meanAmp += curAmp;
       if (peakAmp < curAmp) {
           peakAmp = curAmp;
           peakLagChan = chan;
       }
  }
  ASKAPCHECK(peakLagChan < lags.nelements(), "Empty spectrum is passed to getDelayWithFFT");
  meanAmp -= peakAmp;
  if (lags.nelements() > 1) {
      const double bandwidth = vis.nelements() * itsResolution;
      const double delay =  (static_cast<casacore::Int>(peakLagChan) - static_cast<casacore::Int>(vis.nelements()) / 2)  / bandwidth;
      meanAmp /= double(lags.nelements() - 1);
      ASKAPDEBUGASSERT(meanAmp >= 0.);
      // atan2 is a convenient function to map a ratio of two non-negative numbers to the [0,1] interval 
      itsQuality = atan2(peakAmp, meanAmp) * casacore::C::_2_pi;
      ASKAPDEBUGASSERT(itsQuality >= 0.);
      ASKAPDEBUGASSERT(itsQuality <= 1.);      
      return delay;
  }
  // degenerate case of a single spectral point - unable to estimate delay
  itsQuality = 0.;
  return 0.;
}


/// @brief set new spectral resolution
/// @details The new value will apply to all subsequent calculations
/// @param[in] resolution the spectral resolution in Hz
void DelayEstimator::setResolution(const double resolution) 
{
  itsResolution = resolution;
}

} // namespace scimath

} // namespace askap
