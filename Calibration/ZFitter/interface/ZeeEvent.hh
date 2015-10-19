#ifndef zeeevent_hh
#define zeeevent_hh

/// class ZeeEvent ZeeEvent.h "interface/ZeeEvent.h"

class ZeeEvent{
 public:
  float energy_ele1;
  float energy_ele2;
  float pt1;//This is needed if you want to deal with pt1 + pt2
  float pt2;
  std::vector<float> targetVariable;
  float weight;
  bool  isSwapped;
  float *smearings_ele1, *smearings_ele2;
};

typedef std::vector<ZeeEvent> zee_events_t;
#endif
