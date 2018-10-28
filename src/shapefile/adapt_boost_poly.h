#pragma once

#include "opencv_compat.h"
#include "boost_cv_compat.h"
#include "polygon.h"

/* Boost polygon traits for c_ply */
/* namespace boost {
  template <typename Pt> range_iterator<c_ply<Pt> > {
    typedef vertex_iterator<Pt> type;
  };
} // end namespace boost
*/

namespace boost { namespace geometry {

  namespace traits {

    /*** ply_vertex as point_concept ***/
    template<typename Pt> struct tag<ply_vertex<Pt> *> {
      typedef point_tag type;
    };

    template <typename Pt> struct dimension<ply_vertex<Pt> *>
      : boost::mpl::int_<2> { };

    template <typename Pt> struct access<ply_vertex<Pt> *, 0> {
      static inline typename ply_vertex<Pt>::coordinate_type
        get(ply_vertex<Pt>*p) { return p->getPos().x; }
    };
    template <typename Pt> struct access<ply_vertex<Pt> *, 1> {
      static inline typename ply_vertex<Pt>::coordinate_type
        get(ply_vertex<Pt>*p) { return p->getPos().y; }
    };

    /*** c_ply as ring_concept ***/

    template <typename Pt> struct tag<c_ply<Pt> > {
      typedef ring_tag type;
    };
    template <typename Pt> struct closure<c_ply<Pt> > {
      static const closure_selector value = closed;
    };

    /*** c_polygon as polygon_concept ***/

    template <typename Pt> struct tag<c_polygon<Pt> > {
      typedef polygon_tag type;
    };

    template <typename Pt> struct exterior_ring<c_polygon<Pt> > {
      static inline c_ply<Pt> &get(c_polygon<Pt> &poly)
        { return poly.front(); }
      static inline const c_ply<Pt> &get(const c_polygon<Pt> &poly)
        { return poly.front(); }
    };

    template <typename Pt> struct interior_rings<c_polygon<Pt> > {
      static inline typename c_polygon<Pt>::ring_view
        get(const c_polygon<Pt> &poly) { return poly.inner(); }
    };

  } // end namespace traits

  template <typename Pt> struct ring_type<c_polygon<Pt> > {
    typedef c_ply<Pt> type;
  };

  template <typename Pt> struct interior_type<c_polygon<Pt> > {
    typedef typename c_polygon<Pt>::ring_view type;
  };

} } // end namespace boost::geometry

/* Register ply_vertex with boost as a wrapper for its underlying point.  */
namespace boost { namespace polygon {

  /**** ply_vertex ****/

  template <typename Pt> struct geometry_concept<ply_vertex<Pt>*>
  {
    typedef point_concept type;
  };

  template <typename Pt> struct point_traits<ply_vertex<Pt>*>
  {
    typedef typename ply_vertex<Pt>::coordinate_type coordinate_type;
    static inline coordinate_type get(const ply_vertex<Pt>* pt,
        orientation_2d orient) {
      return (orient == HORIZONTAL) ? pt->getPos().x : pt->getPos().y;
    }
  };

  template <typename Pt> struct point_mutable_traits<ply_vertex<Pt>*>
  {
    typedef typename ply_vertex<Pt>::coordinate_type coordinate_type;
    static inline void set(ply_vertex<Pt> *p, orientation_2d orient,
        coordinate_type value)
    {
      switch(orient.to_int()) {
        case VERTICAL: p->setX(value); break;
        case HORIZONTAL: default: p->setY(value); break;
      }
    }
  };

  /**** c_ply ****/

  template <typename Pt> struct geometry_concept<c_ply<Pt> >
  {
    typedef polygon_concept type;
  };

  template<typename Pt> struct polygon_traits_general<c_ply<Pt> >
  {
    typedef typename c_ply<Pt>::coordinate_type coordinate_type;
    typedef typename c_ply<Pt>::const_iterator iterator_type;
    typedef ply_vertex<Pt> *point_type;

    static inline iterator_type
      begin_points(const c_ply<Pt> &t) { return t.cbegin(); }

    static inline iterator_type
      end_points(const c_ply<Pt> &t) { return t.cend(); }

    static inline unsigned int size(const c_ply<Pt> &ply) {
      return static_cast<unsigned int>(ply.getSize());
    }

    static inline winding_direction winding(const c_ply<Pt> &plt) {
      return unknown_winding;
    }
  };

  /**** c_polygon ****/

  template <typename Pt> struct geometry_concept<c_polygon<Pt> >
  {
    typedef polygon_concept type;
  };

  template<typename Pt> struct polygon_traits_general<c_polygon<Pt> >
  {
    typedef typename c_polygon<Pt>::coordinate_type coordinate_type;
    typedef typename c_polygon<Pt>::vertex_type* point_type;
    typedef typename std::vector<point_type>::const_iterator iterator_type;

    static inline iterator_type begin_points(const c_polygon<Pt> &t)
      { return t.all.cbegin(); }

    static inline iterator_type end_points(const c_polygon<Pt> &t)
      { return t.all.cend(); }

    static inline int size(const c_polygon<Pt> &ply) {
      return ply.getSize();
    }

    static inline winding_direction winding(const c_polygon<Pt> &plt) {
      return unknown_winding;
    }
  };

} } // end namespace boost::polygon
