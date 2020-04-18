#!/usr/bin/env python
# encoding: utf-8
"""
search species and reactions in rmg reactions ans species
"""
from IPython.display import display
from rmgpy.molecule  import molecule

def Findsp(sp,spc1):
    for i in range(len(sp)):
        if spc1 == sp[i].label:
            print sp[i].molecule
            MOL=sp[i]
    return MOL

def FindReactionsInig(ig,spc):
    x1=[]
    for rxn in ig.reactions:
        if spc in [x.name for x in rxn.reactants+rxn.products]:
            x1.append(rxn)
            print [y.name for y in rxn.reactants],"->",[y.name for y in rxn.products]," ",rxn.radicalchange
    return x1

def FindReactionRMG(rmgrxn,spc):
    x1=[]
    for rxn in rmgrxn:
        if spc in [x.label for x in rxn.reactants+rxn.products]:
            x1.append(rxn)
            print rxn.index,[y.label for y in rxn.reactants],"->",[y.label for y in rxn.products]
    return x1

def FindReactionByIndex(rxn,num):
    for react in rxn:
        if react.index==num:
            display(react)
            if hasattr(react, 'library'):
                print react.library
            else:
                print "reaction oblect have no library attribute"
            print react.kinetics,"\n"

def FindAllReactionFromLib(rxn,lib):
    i=1
    for react in rxn:
        if hasattr(react, 'library'):
            if react.library==lib:
                print i,react
                display(react)
                print react.kinetics,"\n"
                i=i+1

def DisplayReactions(rxn_list):
    for react in rxn_list:
        print react
        display(react)
        if hasattr(react, 'library'):
            print react.library
        else:
            print "reaction oblect have no library attribute"
        print react.kinetics,"\n"

def FindbyMW(sp,MW,tol):
    MOL=[]
    for i in range(len(sp)):
        if round(sp[i].molecularWeight.value,tol) == MW:
            print sp[i].molecule
            MOL.append(sp[i])
    if len(MOL) == 0:
        print ("could n't find the species")
    else:
        return MOL

def Findbyinchi(sp,inchi):
    MOL=None
    for i,x in enumerate(sp):
        if str(x.molecule[0].toInChI()) == inchi:
            print x.molecule
            print sp[i].thermo.comment
            MOL=x
    if MOL is None:
        print ("could n't find the species")
    else:
        return MOL

def Findbyisomorph(sp,sp1):
    MOL=None
    res=sp1.generate_resonance_structures()
    for i,x in enumerate(sp):
        for s in res:
            if x.isIsomorphic(s):
                print x.molecule
                print sp[i].thermo.comment
                MOL=x
    if MOL is None:
        print ("could n't find the species")
    else:
        return MOL

def Findbysmiles(sp,smiles):
    sp1=Molecule().fromSMILES(smiles)
    MOL=None
    res=sp1.generate_resonance_structures()
    for i,x in enumerate(sp):
        for s in res:
            if x.isIsomorphic(s):
                print x.molecule
                print sp[i].thermo.comment
                MOL=x
    if MOL is None:
        print ("could n't find the species")
    else:
        return MOL

