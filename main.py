import argparse
import json

def main():
    # Parse information from the command line
    #parser = argparse.ArgumentParser(prog='SCAMPy')
    #parser.add_argument("namelist")
    #parser.add_argument("paramlist")
    #args = parser.parse_args()

    #file_namelist = open(args.namelist).read()
    #namelist = json.loads(file_namelist)
    #del file_namelist

    #file_paramlist = open(args.paramlist).read()
    #paramlist = json.loads(file_paramlist)
    #del file_paramlist

    #main1d(namelist, paramlist)
    main_scampify()

    return

def main1d(namelist, paramlist):
    import Simulation1d
    Simulation = Simulation1d.Simulation1d(namelist, paramlist)
    Simulation.initialize(namelist)
    Simulation.run()
    print('The simulation has completed.')

    return

def main_scampify():
    import scampify
    Simulation = scampify.Scampify1d()
    Simulation.run()
    #Simulation.do_updrafts()
    #Simulation.do_environment()
    #Simulation.do_rain()
    #Simulation.plot_all()
    print('The scampified simulation has completed.')

    return

if __name__ == "__main__":
    main()
