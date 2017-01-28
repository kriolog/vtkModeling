import os
import time


def file_names(template):
    i = 0
    while True:
        while os.path.exists(template%(i,)):
            i += 1
        yield template%(i,)

def toggle(start = False):
    while True:
        start = not start
        yield start


class Timer(object):
    def __init__(self):
        self.running = False
        self.balance = 0
        self.started_at = 0
        
    def start(self):
        if not self.running:
            self.started_at = time.time()
            self.running = True
        
    def stop(self):
        if self.running:
            self.balance += time.time() - self.started_at
            self.running = False
        
    def get_time(self):
        if self.running:
            return self.balance + time.time() - self.started_at
        else:
            return self.balance
        
    def reset(self):
        self.balance = 0
        self.started_at = time.time()
        
    def substract(self, time):
        self.balance -= time
        
    def is_running(self):
        return self.running
        
        
        
        
    
