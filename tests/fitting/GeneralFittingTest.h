/// @file
/// 
/// @brief This file contains general fitting tests
/// @details One of the examples is an equation, which appear in
/// the gain calibration problem (Vmeas= g1*conj(g2)*Vtrue).
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

#ifndef GENERAL_FITTING_TEST_H
#define GENERAL_FITTING_TEST_H


#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/Arrays/Matrix.h>

#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/scimath/fitting/GenericNormalEquations.h>

#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

#include <cppunit/extensions/HelperMacros.h>

#include <cmath>

using std::abs;

namespace askap
{

namespace scimath
{

  class GeneralFittingTest : public CppUnit::TestFixture {

     CPPUNIT_TEST_SUITE(GeneralFittingTest);
     CPPUNIT_TEST(testRealEquationSVD);
     CPPUNIT_TEST(testRealEquationLSQR);
     CPPUNIT_TEST(testComplexEquationSVD);
     CPPUNIT_TEST(testComplexEquationLSQR);
     CPPUNIT_TEST_SUITE_END();

  public:
     void setUp() {
        itsNAnt = 6;
        const casa::uInt nBaselines = itsNAnt*(itsNAnt-1)/2;
        itsBaselines.resize(nBaselines);
        for (casa::uInt ant1=0,baseline=0;ant1<itsNAnt;++ant1) {
             for (casa::uInt ant2=ant1+1;ant2<itsNAnt;++ant2,++baseline) {
                  ASKAPASSERT(baseline<nBaselines);
                  itsBaselines[baseline].first=ant1;
                  itsBaselines[baseline].second=ant2;
             }
        }
        
     }

     void testRealEquationSVD() {
         testRealEquation("SVD");
     }

     void testRealEquationLSQR() {
         testRealEquation("LSQR");
     }
     
     void testComplexEquationSVD() {
         testComplexEquation("SVD");
     }

     void testComplexEquationLSQR() {
         testComplexEquation("LSQR");
     }

  protected:
     /// @brief a helper class to get complex sequence from two real sequences.
     /// @details This helper class acts as an iterator over a complex-valued
     /// sequence. It is initialized with two iterators over two real-valued
     /// sequences, which represent real and imaginary part of the output,
     /// respectively. No boundary checks are done.
     template<typename RealIter, typename ImagIter>
     struct ComplexSequenceIterator {
         /// @brief constructor - remembers iterators
         /// @param[in] ri real part iterator
         /// @param[in] ii imaginary part iterator 
         ComplexSequenceIterator(const RealIter &ri, const ImagIter &ii) :
                  itsRealIter(ri), itsImagIter(ii) {}
         
         /// @brief read-only access method (enough for now)
         /// @return a reference to current element
         casa::Complex operator*() const {
             return casa::Complex(*itsRealIter,*itsImagIter);
         }
         
         /// @brief advance operator
         /// @return reference to itself
         ComplexSequenceIterator& operator++() {
             ++itsRealIter; ++itsImagIter;
             return *this;
         }
     private:
         /// @brief real part iterator
         RealIter itsRealIter;
         /// @brief imaginary part iterator
         ImagIter itsImagIter;
     }; 
     
     /// @brief helper function to create complex sequence iterator
     /// @details This is necessary to avoid specifying template arguments
     /// all the time. These types will be deduced from the method's
     /// arguments types.
     /// @param[in] ri real part iterator
     /// @param[in] ii imaginary part iterator
     template<typename RealIter, typename ImagIter>
     static ComplexSequenceIterator<RealIter, ImagIter> makeComplex(const RealIter &ri,
              const ImagIter &ii) {
        return ComplexSequenceIterator<RealIter,ImagIter>(ri,ii);
     }
  
     /// @brief create parameter object from a sequence of values
     /// @details This method reads itsNAnt elements from the specified
     /// sequence (general iterator semantics) and packs them into 
     /// a Param object. Template is used to allow various types of
     /// input iterators to be used. Because Param::add is an overloaded
     /// method, the template approach also ensures that the correct 
     /// add method will be used for each individual type of the input
     /// values (i.e. complex or real will be determined automatically)
     /// @param[in] iter input iterator of type Iter (template parameter).
     /// The method will read exactly itsNAnt values.
     /// @param[out] params a reference to Params object to fill (will be reset) 
     template<typename Iter> 
     void createParams(const Iter &iter, Params &params) const {
          params.reset();
          Iter cIter(iter);
          for (casa::uInt ant=0; ant<itsNAnt; ++ant,++cIter) {
               params.add(parName(ant),*cIter);
          }
     }
  
     /// @brief obtain parameter name from antenna number
     /// @details
     /// @param[in] ant antenna number
     /// @return string with parameter name
     std::string parName(casa::uInt ant) const {
         ASKAPASSERT(ant<itsNAnt);
         return std::string("gain.")+utility::toString(ant);
     }
     
     /// @brief obtain antenna number from the parameter name
     /// @details Exception is thrown, if parameter name is incompatible
     /// @param[in] par parameter name
     /// @return antenna number [0..itsNAnt-1]
     casa::uInt antNumber(const std::string &par) const {
         ASKAPASSERT(par.find("gain.") == 0 && par.size()>5);
         casa::uInt ant = utility::fromString<casa::uInt>(par.substr(5));
         ASKAPASSERT(ant<itsNAnt);
         return ant;
     }
     
     /// @brief predict products from itsTrueGains
     /// @details This variant implements equation working with
     /// real-valued gains. Predicted values are stored in itsRealMeasuredValues
     void predictReal() {
         itsRealMeasuredValues.resize(itsBaselines.size());
         for (size_t baseline=0; baseline<itsBaselines.size(); ++baseline) {
               itsRealMeasuredValues[baseline] = 
                  itsTrueGains.scalarValue(parName(itsBaselines[baseline].first))*
                  itsTrueGains.scalarValue(parName(itsBaselines[baseline].second));
         }
     }
     
     /// @brief calculate normal equations at itsGuessedGains
     /// @details This variant calculates normal equations for the case
     /// of real-valued gains (i.e. just amplitudes). 
     /// @param[in] ne a reference to normal equations object
     void calcEquationsReal(GenericNormalEquations &ne) {
         ASKAPASSERT(itsBaselines.size() == itsRealMeasuredValues.size());
         casa::Matrix<double> derivatives(itsBaselines.size(),itsNAnt,0.);
         for (size_t baseline=0; baseline<itsBaselines.size(); ++baseline) {
              ASKAPASSERT(itsBaselines[baseline].first<itsNAnt);
              ASKAPASSERT(itsBaselines[baseline].second<itsNAnt);
              derivatives(baseline,itsBaselines[baseline].first) = 
                    itsGuessedGains.scalarValue(parName(itsBaselines[baseline].second));
              derivatives(baseline,itsBaselines[baseline].second) = 
                    itsGuessedGains.scalarValue(parName(itsBaselines[baseline].first));
         }
         DesignMatrix designMatrix; // params: itsGuessedGaussian
         for (casa::uInt ant=0; ant<itsNAnt; ++ant) {
              designMatrix.addDerivative(parName(ant),derivatives.column(ant));
         }
         casa::Vector<double> residual(itsRealMeasuredValues.copy());
         for (size_t baseline=0; baseline<itsBaselines.size(); ++baseline) {
              residual[baseline] -= 
                 itsGuessedGains.scalarValue(parName(itsBaselines[baseline].first))*
                 itsGuessedGains.scalarValue(parName(itsBaselines[baseline].second)); 
         }
         designMatrix.addResidual(residual,casa::Vector<double>(residual.size(),1.));
         ne.add(designMatrix);
     }
       
     /// @brief predict products from itsTrueGains
     /// @details This variant implements equation working with
     /// complex gains. Predicted values are stored in itsComplexMeasuredValues
     void predictComplex() {
         itsComplexMeasuredValues.resize(itsBaselines.size());
         for (size_t baseline=0; baseline<itsBaselines.size(); ++baseline) {
               itsComplexMeasuredValues[baseline] = 
                  itsTrueGains.complexValue(parName(itsBaselines[baseline].first))*
                  conj(itsTrueGains.complexValue(parName(itsBaselines[baseline].second)));
         }
     }

     /// @brief calculate normal equations at itsGuessedGains
     /// @details This variant calculates normal equations for the case
     /// of complex-valued gains (i.e. amplitudes and phases) 
     /// @param[in] ne a reference to normal equations object
     void calcEquationsComplex(GenericNormalEquations &ne) {
         ASKAPASSERT(itsBaselines.size() == itsComplexMeasuredValues.size());
         
         // the first axis has double length because each pair of consequitive
         // elements corresponds to real and imaginary part of the complex-valued
         // gradient. The second axis distinguishes between derivatives by
         // real part and derivatives by imaginary part of the appropriate 
         // gain coefficient.
         casa::Cube<double> derivatives(itsBaselines.size()*2+1,2,itsNAnt,0.);
         casa::Vector<double> residual(itsBaselines.size()*2+1);
         for (size_t baseline=0; baseline<itsBaselines.size(); ++baseline) {
              ASKAPASSERT(itsBaselines[baseline].first<itsNAnt);
              ASKAPASSERT(itsBaselines[baseline].second<itsNAnt);
              casa::Complex g2 = itsGuessedGains.complexValue(parName(
                                        itsBaselines[baseline].second));
              casa::Complex g1 = itsGuessedGains.complexValue(parName(
                                        itsBaselines[baseline].first));
              // d/dRe(g1)=conj(g2)
              derivatives(baseline*2,0,itsBaselines[baseline].first) = real(g2);
              derivatives(baseline*2+1,0,itsBaselines[baseline].first) = -imag(g2);
              // d/dIm(g1)=i*conj(g2)
              derivatives(baseline*2,1,itsBaselines[baseline].first) = imag(g2);
              derivatives(baseline*2+1,1,itsBaselines[baseline].first) = real(g2);
              // d/dRe(g2)=g1
              derivatives(baseline*2,0,itsBaselines[baseline].second) = real(g1);
              derivatives(baseline*2+1,0,itsBaselines[baseline].second) = imag(g1);
              // d/dIm(g2)=-i*g1
              derivatives(baseline*2,1,itsBaselines[baseline].second) = imag(g1);
              derivatives(baseline*2+1,1,itsBaselines[baseline].second) = -real(g1);
 
              const casa::Complex resBuf = itsComplexMeasuredValues[baseline] - 
                                             g1*conj(g2);
              residual[baseline*2] = real(resBuf);
              residual[baseline*2+1] = imag(resBuf);
         }
         /*
         residual[itsBaselines.size()*2] = -imag(itsGuessedGains.complexValue(
                                                 parName(0)));
         derivatives(itsBaselines.size()*2,1,0) = 1.;
         */
         const casa::Complex refGain = itsGuessedGains.complexValue(parName(0));
         residual[itsBaselines.size()*2] = -imag(refGain*itsRefPhase);
         derivatives(itsBaselines.size()*2,0,0) = imag(itsRefPhase);
         derivatives(itsBaselines.size()*2,1,0) = real(itsRefPhase);
         
         DesignMatrix designMatrix; // params: itsGuessedGains;
         for (casa::uInt ant=0; ant<itsNAnt; ++ant) {
              designMatrix.addDerivative(parName(ant),derivatives.xyPlane(ant));
         }
         designMatrix.addResidual(residual,casa::Vector<double>(residual.size(),1.));
         ne.add(designMatrix);
     }

     void testRealEquation(const std::string &solverType) {
         const double trueGains[] = {0.9,1.1,1.2,0.6,1.3,1.05};
         createParams(&trueGains[0],itsTrueGains);
         predictReal();
         const double guessedGains[] = {1.,1.,1.,1.,1.,1.};
         createParams(&guessedGains[0],itsGuessedGains);
         GenericNormalEquations ne;
         calcEquationsReal(ne);
         Quality q;
         LinearSolver solver;
         solver.addNormalEquations(ne);
         solver.setAlgorithm(solverType);
         solver.solveNormalEquations(itsGuessedGains,q);
         for (casa::uInt ant=0;ant<itsNAnt;++ant) {
              CPPUNIT_ASSERT(abs(itsGuessedGains.scalarValue(parName(ant))-
                                 trueGains[ant])<0.05);
         }
     }

     void testComplexEquation(const std::string &solverType) {
         const double trueGainsReal[] = {0.9,1.1,1.2,0.6,1.3,1.05};
         const double trueGainsImag[] = {0.1,-0.1,0.,0.1,-0.1,0.};

         createParams(makeComplex(&trueGainsReal[0], &trueGainsImag[0]),
                      itsTrueGains);
         // correction factor to apply a phase shift to get a correct
         // absolute phase, which we can't determine in the calibration solution.
         // Antenna 0 is a reference.
         itsRefPhase = casa::polar(casa::Float(1.),
                             -casa::arg(itsTrueGains.complexValue(parName(0))));

         predictComplex();
         const double guessedGainsReal[] = {1.,1.,1.,1.,1.,1.};
         const double guessedGainsImag[] = {0.,0.,0.,0.,0.,0.};
         createParams(makeComplex(&guessedGainsReal[0], &guessedGainsImag[0]),
                      itsGuessedGains);
         for (size_t iteration=0; iteration<5; ++iteration) {
              GenericNormalEquations ne;
              calcEquationsComplex(ne);
              Quality q;
              LinearSolver solver;
              solver.addNormalEquations(ne);
              solver.setAlgorithm(solverType);
              solver.solveNormalEquations(itsGuessedGains,q);
         }
         for (casa::uInt ant=0;ant<itsNAnt;++ant) {
             const std::string name = parName(ant);
             CPPUNIT_ASSERT(abs(itsTrueGains.complexValue(name)-
                      itsGuessedGains.complexValue(name))<1e-6);
         }
     }

  private:
     /// @brief guessed parameters
     Params itsGuessedGains;
     /// @brief true parameters
     Params itsTrueGains;
     /// @brief reference phase
     /// @detail A complex number with amplitude of 1 and the phase
     /// equivalent to the phase of the gain of a reference antenna
     /// (we can't determine absolute phase and have to adopt something)
     casa::Complex itsRefPhase;
     /// @brief data sample playing a role of measured real data
     casa::Vector<double> itsRealMeasuredValues;
     /// @brief data sample playing a role of measured complex data
     casa::Vector<casa::Complex> itsComplexMeasuredValues;
     /// @brief first and second "antenna" corresponding to a baseline
     std::vector<std::pair<casa::uInt, casa::uInt> > itsBaselines;
     /// @brief number of antennae
     casa::uInt itsNAnt;
};

} // namespace scimath

} // namespace askap

#endif // #ifndef GENERAL_FITTING_TEST_H
