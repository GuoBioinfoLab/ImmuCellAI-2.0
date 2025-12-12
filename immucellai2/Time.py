#!/usr/bin/python3
import time
import datetime


class CLASS_FOR_TIME(object):
   def __init__(self, recordname="hahaha"):
      self.NowTime = datetime.datetime.now
      self.StartTime = self.NowTime()
      self.recordname = recordname
      self.ShowTime()

   def ShowTime(self):
      print("[%s] Curent time: "%self.recordname, self.NowTime())

   def ShowCostTime(self):
      # str(CostTime)
      self.ShowTime()
      self.EndTime = self.NowTime()
      CostTime = self.EndTime - self.StartTime
      print("[%s] Cost time: "%self.recordname, CostTime)
       
