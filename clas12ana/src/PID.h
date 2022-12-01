#ifndef PID_h
#define PID_h

#include "SIDISParticle.h"
#include "SIDISParticlev1.h"
using namespace std;

typedef std::map<int, SIDISParticle*> type_map_part;

class PID{
 public:
  PID();
  
  type_map_part applyCuts(type_map_part&);

  int countPID(type_map_part&, int);
};

#endif
