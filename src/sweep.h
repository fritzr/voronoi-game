#pragma once

#include <queue>

template<class Event, class EventCompare, class EventId=int>
struct event_traits
{
  typedef Event   event_type;
  typedef EventId event_id_type;

  inline static event_id_type const& get_type(event_type const& e)
  { return e.type; }
};

template<
  class Value, class Event, class Solution,
  class EventTraits=event_traits<Event, std::less<Event>, int>,
  class EventContainer=std::vector<Event>,
  class SolutionContainer=std::vector<Solution>
  >
struct sla_traits
{
  typedef Value value_type;
  typedef EventTraits event_traits;
  typedef typename event_traits::event_type    event_type;
  typedef typename event_traits::event_compare event_compare;
  typedef typename event_traits::event_id_type event_id_type;
  typedef EventContainer                       event_container;
  typedef std::priority_queue<event_type, event_container, event_compare>
    queue_type;

  typedef Solution solution_type;
  typedef SolutionContainer solution_container;
};

template<class Traits>
class SweepLineAlgorithm
{
public:
  // Typedefs
  typedef Traits traits;
  typedef typename traits::value_type      value_type;
  typedef typename traits::event_type      event_type;
  typedef typename traits::event_compare   event_compare;
  typedef typename traits::event_id_type   event_id_type;
  typedef typename traits::event_container event_container;
  typedef typename traits::queue_type      queue_type;

  typedef typename traits::solution_type      solution_type;
  typedef typename traits::solution_container solution_container;

  inline static event_id_type const& get_type(event_type const& e) {
    return traits::event_traits::get_type(e);
  }

protected:
  // Members
  queue_type queue_;
  inline queue_type& queue(void) { return queue_; }
  inline queue_type const& queue(void) const { return queue_; }

  solution_container solutions_;
  inline solution_container& solutions(void) { return solutions_; }
  inline solution_container const& solutions(void) const { return solutions_; }

  virtual void handle_event(event_id_type const& ty, event_type const& event)=0;
  virtual void initialize(void) {}
  virtual void finalize(void) {}

public:
  // Methods
  SweepLineAlgorithm(void) {}
  virtual ~SweepLineAlgorithm(void) {}

  virtual void add_event(value_type const& v)=0;

  template<class ValueIter>
    inline void insert(ValueIter begin, ValueIter end) {
      while (begin != end)
        add_event(*begin++);
    }

  void compute(void)
  {
    initialize();
    while (!queue().empty())
    {
      typename queue_type::const_reference event = queue().top();
      handle_event(get_type(event), event);
      queue().pop();
    }
    finalize();
  }
};
