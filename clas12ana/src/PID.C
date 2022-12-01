#include "PID.h"

using namespace std;


PID::PID(){}

type_map_part PID::applyCuts(type_map_part& map){


  type_map_part retMap;

  // Find the scattered electron
  TLorentzVector e;
  for(type_map_part::iterator it = map.begin(); it!=map.end(); ++it){
    int pid = (it->second)->get_property_int(SIDISParticle::part_pid);
    if(pid==11){
      float px = (it->second)->get_property_float(SIDISParticle::part_px);
      float py = (it->second)->get_property_float(SIDISParticle::part_py);
      float pz = (it->second)->get_property_float(SIDISParticle::part_pz);
      float E = (it->second)->get_property_float(SIDISParticle::part_E);
      
      e.SetPxPyPzE(px,py,pz,E);
      break; // Assume only one scattered electron
    }
    else
      continue;
  }

  // Loop over particles within the particle map
  for(type_map_part::iterator it = map.begin(); it!=map.end(); ++it){
    

    int pid = (it->second)->get_property_int(SIDISParticle::part_pid);
    float E = (it->second)->get_property_float(SIDISParticle::part_E);
    float px = (it->second)->get_property_float(SIDISParticle::part_px);  
    float py = (it->second)->get_property_float(SIDISParticle::part_py);  
    float pz = (it->second)->get_property_float(SIDISParticle::part_pz);  
    float chi2 = (it->second)->get_property_float(SIDISParticle::part_chi2);
    float beta = (it->second)->get_property_float(SIDISParticle::part_beta);
    float vz = (it->second)->get_property_float(SIDISParticle::part_vz);

    TLorentzVector part(px,py,pz,E);
    float angle_e = part.Angle(e.Vect())*180/3.14159265;

    // Toss photons near electron (radiative)
    if(pid == 22 && angle_e<8)
      continue;
    
    // Toss low energy photons
    if(pid == 22 && E<0.2)
      continue;

    // Toss bad electron vertex
    if(pid == 11 && (vz<-8 || vz>3))
      continue;

    // Toss bad photon beta
    if(pid == 22 && (beta<0.9 || beta>1.1))
      continue;

    // Insert particle into the retMap
    retMap.insert ( make_pair( (it->first) , (it->second)) );

  }
  
  return retMap;
}


int PID::countPID(type_map_part& map, int pid){
  int count = 0;
  // Loop over particles within the particle map
  for(type_map_part::iterator it = map.begin(); it!=map.end(); ++it){
    if( (it->second)->get_property_int(SIDISParticle::part_pid) == pid)
      count++;
  }
  return count;
}
