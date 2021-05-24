#ifndef GEMOBJECTS_GEMDIGISIMLINK_H
#define GEMOBJECTS_GEMDIGISIMLINK_H

#include <map>
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"

class GEMDigiSimLink {
public:
  GEMDigiSimLink(std::pair<unsigned int, int> digi,
                 Local3DPoint entryPoint,
                 LocalVector momentumAtEntry,
                 float timeOfFlight,
                 float energyLoss,
                 int particleType,
                 unsigned int detUnitId,
                 unsigned int trackId,
                 EncodedEventId eventId,
                 unsigned short processType) {
    _entryPoint = entryPoint;
    _momentumAtEntry = momentumAtEntry;
    _timeOfFlight = timeOfFlight;
    _energyLoss = energyLoss;
    _particleType = particleType;
    _detUnitId = detUnitId;
    _trackId = trackId;
    _eventId = eventId;
    _processType = processType;
    _digi = digi;
  }

  GEMDigiSimLink() { ; }

  ~GEMDigiSimLink() { ; }

  unsigned int getStrip() const { return _digi.first; }
  int getBx() const { return _digi.second; }
  Local3DPoint getEntryPoint() const { return _entryPoint; }
  LocalVector getMomentumAtEntry() const { return _momentumAtEntry; }
  float getTimeOfFlight() const { return _timeOfFlight; }
  float getEnergyLoss() const { return _energyLoss; }
  int getParticleType() const { return _particleType; }
  unsigned int getDetUnitId() const { return _detUnitId; }
  unsigned int getTrackId() const { return _trackId; }
  EncodedEventId getEventId() const { return _eventId; }
  unsigned short getProcessType() const { return _processType; }

  inline bool operator<(const GEMDigiSimLink& other) const { return getStrip() < other.getStrip(); }

private:
  std::pair<unsigned int, int> _digi;

  Local3DPoint _entryPoint;
  LocalVector _momentumAtEntry;
  float _timeOfFlight;
  float _energyLoss;
  int _particleType;
  unsigned int _detUnitId;
  unsigned int _trackId;
  EncodedEventId _eventId;
  unsigned short _processType;
};
#endif