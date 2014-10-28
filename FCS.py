# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 18:34:54 2014
Creates an FCS object given a filename

@author: David Ng, MD
"""
"""Installed Packages"""
import numpy as np
import pandas as pd
"""Built in packages"""
from re import compile, findall
from datetime import datetime
from warnings import warn
from struct import calcsize, unpack
from os.path import basename, relpath, dirname

from FCS_subroutines.loadFCS import loadFCS
from FCS_subroutines.Process_FCS_Data import Process_FCS_Data
from FCS_subroutines.empty_FCS import empty_FCS
from FCS_subroutines.FCSmeta_to_database import FCSmeta_to_database

class FCS(object):
    """
    This class represents FCS data (Tube+Case information)
    See loadFCS for attribute details
    """
    def __init__(self, version, filepath=None, db=None, **kwargs):
        if filepath is not None and db is not None:
            raise "Must import data from file or db, not both!"
        if filepath is not None:
            self.load_from_file(filepath, version, **kwargs)
        elif db is not None:
            self.load_from_db(db)

    def load_from_file(self, filepath, version, **kwargs):
        """ 
        Import FCS data from file at <filepath>
        nota bene: import_dataframe needs to be explicitly defined for
        data to be loaded into FCS object        
        """
        loadFCS(FCS=self, filepath=filepath, version=version, **kwargs)

    def load_from_db(self, db):
        """ Import FCS data from db <db> """
        raise "Not implemented"

    def meta_to_db(self, db, dir=None, add_lists=False):
        """ Export meta data from FCS object to db """
        FCSmeta_to_database(FCS=self, db=db, dir=dir, add_lists=add_lists)

    def comp_scale_FCS_data(self,compensation_file,
                            saturation_upper_range=1000,
                            rescale_lim=(-0.15, 1), 
                            limits=False, 
                            strict=True, 
                            **kwargs):
        """calls Process_FCS_Data on self (i.e. an FCS object)"""
        Process_FCS_Data(FCS = self, compensation_file = compensation_file, 
                         saturation_upper_range = saturation_upper_range,
                         rescale_lim = rescale_lim,
                         limits = limits,
                         strict = strict,
                         **kwargs)



if __name__ == '__main__':
    filename = "/home/ngdavid/Desktop/MDS_Plates/Hodgkin_Cases_2008_2013/10-06255/10-06255_Hodgkins.fcs"
    filename = "/home/ngdavid/Desktop/MDS_Plates/12-02814/Plate 3/12-02814_C11_C11.fcs"
    temp = FCS(version = '1')
    temp.load_from_file(filename, version=1, import_dataframe=True)
