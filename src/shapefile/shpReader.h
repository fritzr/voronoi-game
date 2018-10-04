//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _SHPREADER_H_
#define _SHPREADER_H_

#include "polygon.h"

//
// this reader reads ESRI shapefile format and convert
// polygons to poly format
//

class ShapeReader 
{
public:

    //read from file
    bool read(const string& filename);

    //get poly read from file
    vector<c_polygon>& getPolyList(){ return m_ply_list; }

    //set tags for getting information from database
    void setBaseElevationTag(const string& tag) { base_elevation_tag=tag; }
    void setHeightTag(const string& tag) { height_tag=tag; }

private:

    void add_ply(c_ply& plys);

    vector<c_polygon> m_ply_list; //a list of polygons
    string base_elevation_tag;
    string height_tag;
};

#endif //_SHPREADER_H_
