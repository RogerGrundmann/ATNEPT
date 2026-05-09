#!/usr/bin/env python

from pyatnept import Neptune

class Model(object):
    """
    ATNEPT Neptune' cell structure in the atmosphere
    """
    def __init__(self, cfg_xml):
        self.nept = Neptune()
        self.config_xml = cfg_xml

    def print_config_nept(self, config_xml):
        print("\n\n\n   Neptune' cell structure in the atmosphere") 
        print("\n\n\n   Neptune atmosphere configuration file name = ", self.config_xml)
        print("   output path is           ", self.nept.output_path.decode('utf-8'))

    def run_Model_nept(self):
        print("\n   run_Model for the Neptune-Atmosphere code prepared")
        self.nept.run()
        print("\n    successfully terminated Neptune-Atmosphere code")
        print("\n")

nept = Model("config_atnept.xml")
nept.print_config_nept("config_atnept.xml")
nept.run_Model_nept()
