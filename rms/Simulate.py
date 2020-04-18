#!/usr/bin/env python
# encoding: utf-8
"""
Simulate rms
"""
def rmsconttp(Temp,Patm,Time):
    Ppa=Patm*101325
    initialconds = {"T":Temp,"P":Ppa,"N2":0.78504,"ttDMP":0.018692,"O2":0.19626}
    domain,y0 = rms.ConstantTPDomain(phase=ig,initialconds=initialconds)
    react = rms.Reactor(domain,y0,(0.0,Time))
    sol = de.solve(react.ode,de.CVODE_BDF(),abstol=1e-20,reltol=1e-3)
    #reltol=1e-3 to make run faster eltol=1e-8 is the default
    sim = rms.Simulation(sol,domain)
    return sim

def rmstp(Temp,Patm):
    Ppa=Patm*101325
    initialconds = {"T":Temp,"P":Ppa,"N2":0.78504,"ttDMP":0.018692,"O2":0.19626}
    domain,y0 = rms.ConstantVDomain(phase=ig,initialconds=initialconds)
    react = rms.Reactor(domain,y0,(0.0,0.005))
    sol = de.solve(react.ode,de.CVODE_BDF(),abstol=1e-20,reltol=1e-8)
    #reltol=1e-3 to make run faster reltol=1e-8 is the default
    sim = rms.Simulation(sol,domain)
    return sim