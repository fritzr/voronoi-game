//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "shpReader.h"
#include "shapelib/shapefil.h"

#include <cassert>
#include <iostream>
using namespace std;

typedef unsigned int uint;

//
// read from files (both ESRI shp and db files)
//
bool ShapeReader::read(const string& filename)
{
  SHPHandle  hSHP;
  hSHP = SHPOpen( filename.c_str(), "rb" );

  //something is wrong
  if( hSHP == NULL )
  {
    cout<<"- Error: Unable to open: "<<filename<<endl;
    return false;
  }

  //Get info
  int    nShapeType, nEntities;
  double   adfMinBound[4], adfMaxBound[4];
  SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

  //reserve space
  m_ply_list.reserve(nEntities);

  //  Skim over the list of shapes, printing all the vertices
  for(int i=0; i<nEntities; i++ )
  {

    //
    SHPObject * psShape = SHPReadObject( hSHP, i );

    //make sure psShape is not null
    assert(psShape);

    //some format check..
    if( psShape->nParts > 0 && psShape->panPartStart[0] != 0 )
    {
      cerr<<"panPartStart[0] = "<<psShape->panPartStart[0]
        <<", not zero as expected.\n";
    }

    //read each part
    int iPart = 1;

    //
    c_ply ply(c_ply::POUT);
    ply.beginPoly();
    ply_extra_info& extra=ply.extra();
    extra.box[0]=psShape->dfXMin;
    extra.box[1]=psShape->dfXMax;
    extra.box[2]=psShape->dfYMin;
    extra.box[3]=psShape->dfYMax;
    extra.has_box=true;

    for( int j = 0; j < psShape->nVertices; j++ )
    {
      string pszPartType;

      if( j == 0 && psShape->nParts > 0 )
        pszPartType = SHPPartTypeName( psShape->panPartType[0] );

      //this defines a new loop...
      if( iPart < psShape->nParts && psShape->panPartStart[iPart] == j )
      {
        pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
        iPart++;
        ply.endPoly(true); //remove duplicate
        add_ply(ply);

        //create a new ply
        ply=c_ply(c_ply::PIN);
        ply.beginPoly();
      }

      //Points may duplicate....
      ply.addVertex(psShape->padfX[j],psShape->padfY[j], true);
    }

    ply.endPoly(true); //use true to set the flag to remove duplicate vertices at the end
    add_ply(ply);

    SHPDestroyObject( psShape );
  }

  SHPClose( hSHP );

  //done reading shp file

  //
  //reading from database now.
  //

  DBFHandle hDBF = DBFOpen( filename.c_str(), "rb" );

  //check
  if( hDBF == NULL )
  {
    cerr<<"! Error: DBFOpen("<< filename<<") failed."<<endl;
    return false;
  }

  //check
  if( DBFGetFieldCount(hDBF) == 0 )
  {
    cerr<<"! Warning: There are no fields in this table!"<<endl;
    return true;
  }

  //check more
  if((uint)DBFGetRecordCount(hDBF)!=m_ply_list.size()){
    cerr<<"! Error: Number of rows ("<<DBFGetRecordCount(hDBF)
      <<") in the database does not match to the number of buildings ("
      <<m_ply_list.size()<<")"<<endl;
  }

  //used for reading from db
  int nWidth, nDecimals;

  //position of these fields in database...
  uint base_ele_pos=UINT_MAX;
  uint heigh_pos=UINT_MAX;
  for(int i = 0; i < DBFGetFieldCount(hDBF); i++ )
  {
    char  szTitle[12];
    DBFFieldType eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
    if(eType!=FTInteger && eType!=FTDouble) continue; //not good...
    if(base_elevation_tag==szTitle) base_ele_pos=i;
    if(height_tag==szTitle) heigh_pos=i;
  }

  //check if we can get the field ids for the information we wanted
  if(base_ele_pos==UINT_MAX)
    cerr<<"! Warning: no base elevation found in the database"<<endl;

  if(heigh_pos==UINT_MAX)
    cerr<<"! Warning: no height found in the database"<<endl;

  //
  //read data from each row
  //
  for( int row = 0; row<DBFGetRecordCount(hDBF); row++ )
  {

    ply_extra_info& extra=m_ply_list[row].front().extra();

    if(base_ele_pos!=UINT_MAX)
    {
      //
      bool is_null=DBFIsAttributeNULL( hDBF, row, base_ele_pos );
      if(is_null){
        cerr<<"! Warning: Database for ("<<base_elevation_tag<<") is null at row: "<<row<<endl;
        continue;
      }
      extra.base_height=DBFReadDoubleAttribute( hDBF, row, base_ele_pos );
    }

    if(heigh_pos!=UINT_MAX)
    {
      //
      bool is_null=DBFIsAttributeNULL( hDBF, row, heigh_pos );
      if(is_null){
        cerr<<"! Warning: Database for ("<<height_tag<<") is null at row: "<<row<<endl;
        continue;
      }

      extra.height=DBFReadDoubleAttribute( hDBF, row, heigh_pos );
    }

  }//end row

  DBFClose( hDBF );

  return true;
}


void ShapeReader::add_ply(c_ply& ply)
{
  if(ply.getType()==c_ply::POUT){
    //create a new polygon
    c_polygon polygon;
    polygon.push_back(ply);
    m_ply_list.push_back(polygon);
  }
  else{ //c_ply::PIN
    //add to the last polygon
    //make sure that this hole is inside the last polygon...
    c_polygon& polygon=m_ply_list.back();
    Point2d pt=ply.findEnclosedPt();
    bool b=polygon.enclosed(pt);
    if(b)
      polygon.push_back(ply);
    else
      cerr<<"! WARNING: Input Error: hole is ignored due improper nesting."<<endl;
  }
}


