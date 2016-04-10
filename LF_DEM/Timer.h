/**
 \class Simulation
 \brief Class launching the simulation by setting up the System class and performing predefined shear evolution
 \author Ryohei Seto
 \author Romain Mari
 */
//
//  Timer.h
//  LF_DEM
//
//  Created by Ryohei Seto and Romain Mari on 04/09/2016.
//  Copyright (c) 2016 Ryohei Seto and Romain Mari. All rights reserved.
//
#ifndef __LF_DEM__Timer__
#define __LF_DEM__Timer__
#include <memory>
#include <string>
#include <vector>
#include <stdexcept>

class Clock
{
protected:
  double next_time;
  double time_step;
  std::string _label;
  bool strain;
public:
  std::string label() {
    return _label;
  }
  double nextTime() {
    return next_time;
  }
  bool is_strain(){
    return strain;
  }
  virtual void tick() {};
};

class LinearClock : public Clock
{
private:

public:
  LinearClock(double stop, double time_step, std::string label)
  {
    next_time = 0;
    time_step = time_step;
    _label = label;
  }
  virtual void tick() {
    next_time = exp(log(next_time)+time_step);
  }
};

class LogClock : public Clock
{
private:

public:
  LogClock(double start, double stop, double nb_step, std::string label)
  {
    next_time = start;
    time_step = (log(stop) - log(start))/nb_step;
    _label = label;
  }
  virtual void tick() {
    next_time = exp(log(next_time)+time_step);
  }
};


class TimeKeeper
{
private:
  std::vector <std::unique_ptr<Clock>> clocks;
public:
  void addClock(LinearClock c){
    clocks.push_back(std::unique_ptr<Clock>(new LinearClock(c)));
  }
  void addClock(LogClock c){
    clocks.push_back(std::unique_ptr<Clock>(new LogClock(c)));
  }

  std::pair<double,std::string> nextTime(double current_time) {
    if (clocks.size()==0) {
      throw std::runtime_error( " TimeKeeper::nextStrain() : No clocks! ");
    }
    double next_t = current_time;
    std::string next_name = "";
    for (auto &c : clocks) {
      if (!c->is_strain() && c->nextTime() < next_t) {
        next_t = c->nextTime();
        next_name = c->label();
      }
    }
    return std::make_pair(next_t,next_name);
  }

  std::pair<double,std::string> nextStrain(double current_strain) {
    if (clocks.size()==0) {
      throw std::runtime_error( " TimeKeeper::nextStrain() : No clocks! ");
    }
    if (current_strain<0) {
      throw std::runtime_error( " TimeKeeper::nextStrain() : current_strain must be >0");
    }
    double next_t = current_strain;
    std::string next_name = "";
    for (auto &c : clocks) {
      if (c->is_strain() && c->nextTime() < next_t) {
        next_t = c->nextTime();
        next_name = c->label();
      }
    }
    return std::make_pair(next_t,next_name);
  }

  void updateClocks(double time) {
    for (auto &c : clocks) {
      if (c->nextTime() <= time) {
        c->tick();
      }
    }
  }

  // TimeKeeper();
  // ~TimeKeeper();
};

#endif /* defined(__LF_DEM__Timer__) */
