/// @file
///
/// Holds the parameters in a CASA table
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///
#include <askap/scimath/fitting/ParamsCasaTable.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/Axes.h>

#include <std::string>
#include <iostream>

#include <casacore/casa/aips.h>
#include <casacore/tables/Tables/TableLocker.h>
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/ColumnDesc.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/TaQL/ExprNode.h>
#include <casacore/tables/TaQL/ExprNodeSet.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/ArrayUtil.h>

#include <casacore/casa/Arrays/Slice.h>


#include <casacore/casa/Utilities/Regex.h>
#include <casacore/casa/Utilities/GenSort.h>
#include <casacore/casa/BasicMath/Math.h>

#include <askap/askap/AskapError.h>

using namespace askap;

using namespace casa;

namespace askap
{
  namespace scimath
  {
    /// The column identifiers
    /// @name  Column Identifiers
    //@{
    const String colName("NAME");
    const String colValues("VALUES");
    const String colAxes("AXES");
    const String colStart("AXESSTART");
    const String colEnd("AXESEND");
    const String colDomain("DOMAIN");
    const String colDomainStart("DOMAINSTART");
    const String colDomainEnd("DOMAINEND");
    const String colFree("FREE");
    //@}

    ParamsCasaTable::ParamsCasaTable(const std::string& tablename, bool exists)
    {
      if(exists)
      {
        openTable(tablename);
      }
      else
      {
        createTable(tablename);
      }
    }

    void ParamsCasaTable::createTable(const std::string& tablename)
    {
      itsTableName=tablename;
      itsTableDesc.addColumn (ScalarColumnDesc<String>(colName));
      itsTableDesc.addColumn (ArrayColumnDesc<String>(colAxes));
      itsTableDesc.addColumn (ArrayColumnDesc<double>(colStart,1));
      itsTableDesc.addColumn (ArrayColumnDesc<double>(colEnd,1));
      itsTableDesc.addColumn (ArrayColumnDesc<String>(colDomain));
      itsTableDesc.addColumn (ArrayColumnDesc<double>(colDomainStart,1));
      itsTableDesc.addColumn (ArrayColumnDesc<double>(colDomainEnd,1));
      itsTableDesc.addColumn (ArrayColumnDesc<double>(colValues,-1));
      itsTableDesc.addColumn (ScalarColumnDesc<bool>(colFree));

      SetupNewTable newtab(itsTableName, itsTableDesc, Table::New);
      itsTable=Table(newtab,0,False,Table::LocalEndian);
      std::cout << "Successfully created new parameters table " << itsTableName << std::endl;
    }

    void ParamsCasaTable::openTable(const std::string& tablename)
    {
      ASKAPCHECK(Table::isReadable(tablename), "Parameters table " << tablename << " is not readable");

      itsTableName=tablename;
      itsTable=Table(itsTableName);
      std::cout << "Successfully opened existing parameters table " << itsTableName << std::endl;
    }

    ParamsCasaTable::~ParamsCasaTable()
    {
      itsTable.flush(true);
    }

    void ParamsCasaTable::getParameters(Params& ip) const
    {
      Domain null;
      null.add("NULL", 0.0, 0.0);
      getParameters(ip, null);
    }

    void ParamsCasaTable::getParameters(Params& ip, const Domain&) const
    {
      ASKAPCHECK(Table::isReadable(itsTable.tableName()), "Parameters table " << itsTable.tableName() << " is not readable");

      ROScalarColumn<String> nameCol (itsTable, colName);
      ROArrayColumn<double> valCol (itsTable, colValues);
      ROArrayColumn<String> axesCol (itsTable, colAxes);
      ROArrayColumn<double> startCol (itsTable, colStart);
      ROArrayColumn<double> endCol (itsTable, colEnd);
      ROArrayColumn<String> domainCol (itsTable, colDomain);
      ROArrayColumn<double> domainStartCol (itsTable, colDomainStart);
      ROArrayColumn<double> domainEndCol (itsTable, colDomainEnd);
      ROScalarColumn<bool> freeCol (itsTable, colFree);

      ASKAPCHECK(itsTable.nrow()>0, "Parameters table " << itsTable.tableName() << " is empty");

      for (size_t rownr=0; rownr<itsTable.nrow(); ++rownr)
      {
        casacore::String name;
        nameCol.get(rownr, name);
        casacore::Array<double> value;
        valCol.get(rownr, value);

        casacore::Vector<String> axesNames;
        axesCol.get(rownr, axesNames);
        casacore::Vector<double> start;
        startCol.get(rownr, start);
        casacore::Vector<double> end;
        endCol.get(rownr, end);
        Axes ax;
        for (size_t i=0; i<axesNames.nelements(); ++i)
        {
          ax.add(axesNames(i), start(i), end(i));
        }
        ip.add(name, value, ax);
        bool free;
        freeCol.get(rownr, free);
        if(free) 
        {
          ip.free(name);
        }
        else 
        {
          ip.fix(name);
        }
      }
    };

    casacore::Vector<casacore::String> ParamsCasaTable::toCasaString(const std::vector<std::string>& s)
    {
      casacore::Vector<casacore::String> result(s.size());
      for (size_t i=0; i<s.size(); ++i)
      {
        result(i)=s[i];
      }
      return result;
    }

    std::vector<std::string> ParamsCasaTable::toStdString(const casacore::Vector<casacore::String>& s)
    {
      std::vector<std::string> result(s.nelements());
      for (size_t i=0; i<s.nelements(); ++i)
      {
        result[i]=s(i);
      }
      return result;
    }

    void ParamsCasaTable::setParameters(const Params& ip)
    {
      Domain null;
      null.add("NULL", 0.0, 0.0);
      setParameters(ip, null);
    }

    void ParamsCasaTable::setParameters(const Params& ip, const Domain& domain)
    {
      itsTable.reopenRW();
      TableLocker locker(itsTable, FileLocker::Write);
      ScalarColumn<String> nameCol (itsTable, colName);
      ArrayColumn<double> valCol (itsTable, colValues);
      ArrayColumn<String> axesCol (itsTable, colAxes);
      ArrayColumn<double> startCol (itsTable, colStart);
      ArrayColumn<double> endCol (itsTable, colEnd);
      ArrayColumn<String> domainCol (itsTable, colDomain);
      ArrayColumn<double> domainStartCol (itsTable, colDomainStart);
      ArrayColumn<double> domainEndCol (itsTable, colDomainEnd);
      ScalarColumn<bool> freeCol (itsTable, colFree);

      std::vector<std::string> names(ip.names());
      int rownr=itsTable.nrow();

      for (std::vector<std::string>::iterator it=names.begin();it!=names.end();it++)
      {
        itsTable.addRow();
        nameCol.put(rownr, *it);
        valCol.put(rownr, ip.value(*it));

        Axes ax(ip.axes(*it));
        axesCol.put(rownr, toCasaString(ax.names()));
        startCol.put(rownr, casacore::Vector<double>(ax.start()));
        endCol.put(rownr, casacore::Vector<double>(ax.end()));

        domainCol.put(rownr, toCasaString(domain.names()));
        domainStartCol.put(rownr, casacore::Vector<double>(domain.start()));
        domainEndCol.put(rownr, casacore::Vector<double>(domain.end()));

        freeCol.put(rownr, ip.isFree(*it));

        rownr++;
      }

    }

  }
}
