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
#include <map>
#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>

class Clock
{
protected:
	double next_time;
	double time_step;
	bool _strain;
public:
	Clock(bool strain):
	_strain(strain) {};

	double nextTime()
	{
		return next_time;
	}
	bool is_strain()
	{
		return _strain;
	}
	virtual void tick() {};
};

class LinearClock : public Clock
{
private:

public:
	LinearClock(double step, bool strain_units, double start=0)
	: Clock(strain_units)
	{
		next_time = start;
		time_step = step;
	}

	virtual void tick() {
		next_time += time_step;
	}
};

class LogClock : public Clock
{
private:

public:
	LogClock(double start, double stop, double nb_step, bool strain_units)
	: Clock(strain_units)
	{
		next_time = start;
		time_step = (log(stop) - log(start))/nb_step;
	}

	virtual void tick()
	{
		next_time = exp(log(next_time)+time_step);
	}
};


class TimeKeeper
{
private:
	std::map <std::string, Clock*> clocks;
public:
	void addClock(std::string label, LinearClock c)
	{
		clocks[label] = new LinearClock(c);
	}

	void addClock(std::string label, LogClock c)
	{
		clocks[label] = new LogClock(c);
	}

	void removeClock()
	{
		clocks.clear();
	}

	std::pair<double,std::string> nextTime() const
	{
		if (clocks.size() == 0) {
			throw std::runtime_error( " TimeKeeper::nextTime() : No clocks! ");
		}
		double next_time = -1;
		std::string next_name = "";
		for (const auto &c : clocks) {
			const auto &label = c.first;
			const auto &clk = c.second;
			if (!clk->is_strain() && (clk->nextTime() < next_time || next_time < 0)) {
				next_time = clk->nextTime();
				next_name = label;
			}
		}
		return std::make_pair(next_time, next_name);
	}

	std::pair<double,std::string> nextStrain() const
	{
		if (clocks.size() == 0) {
			throw std::runtime_error( " TimeKeeper::nextStrain() : No clocks! ");
		}
		double next_time = -1;
		std::string next_name = "";
		for (const auto &c : clocks) {
			const auto &label = c.first;
			const auto &clk = c.second;
			if (clk->is_strain() && (clk->nextTime() < next_time || next_time < 0)) {
				next_time = clk->nextTime();
				next_name = label;
			}
		}
		return std::make_pair(next_time, next_name);
	}

	std::set<std::string> getElapsedClocks(double time_, double strain)
	{
		std::set<std::string> elapsed_clocks;
		for (auto &c : clocks) {
			auto &clk = c.second;
			if (clk->is_strain()) {
				if (clk->nextTime() <= strain+1e-8) {
					clk->tick();
					elapsed_clocks.insert(c.first);
				}
			} else {
				if (clk->nextTime() <= time_+1e-8) {
					clk->tick();
					elapsed_clocks.insert(c.first);
				}
			}
		}
		return elapsed_clocks;
	}
	// TimeKeeper();
	// ~TimeKeeper();
};

#endif /* defined(__LF_DEM__Timer__) */
