

class Scheduler(object):


    def __init__(self):
        
        self.now = self.getCurrentTime()

        self.alarmIndex = []
        self.alarms = {}

        return


    def alarm(self, interval, callback):
        """Call the given callback after the specified time interval
        elapses."""

        from pyre.units.time import second
        alarmTime = self.now + interval/second

        newAlarm = self.Alarm(alarmTime)
        alarm = self.alarms.setdefault(alarmTime, newAlarm)
        alarm.append(callback)

        if alarm is newAlarm:
            self.alarmIndex.append(alarmTime)
            self.alarmIndex.sort(reverse = True)
        
        return


    def poll(self):
        """Call the callbacks for any alarms that have gone off.
        Answer the number of seconds we can sleep until the next
        alarm.  If there are no more alarms, answer None."""

        self.updateInternalClock()

        activeAlarm = self.activeAlarm
        
        if activeAlarm is None:
            return None # sleep indefinitely
        
        while activeAlarm.time <= self.now:
            for callback in activeAlarm:
                callback()
            # activate the next alarm
            activeAlarm = self.popActiveAlarm()
            if activeAlarm is None:
                return None # sleep indefinitely

        return activeAlarm.time - self.now


    # private

    def updateInternalClock(self):
        """Advance our internal clock to the current system time."""
        now = self.getCurrentTime()
        if now < self.now:
            self.clockSetBack(self.now - now)
        self.now = now
        return


    def clockSetBack(self, delta):
        """The system clock was set back; update our internal data
        structures."""
        if not self.alarms:
            return # nothing to do
        self.alarmIndex = []
        for alarm in self.alarms:
            alarm.time -= delta
            self.alarmIndex.append(alarm.time)
        self.alarmIndex.sort(reverse = True)
        return


    def getActiveAlarm(self):
        if self.alarmIndex:
            return self.alarms[self.alarmIndex[-1]]
        return None
    activeAlarm = property(getActiveAlarm)


    def popActiveAlarm(self):
        """Discard the currently active alarm.  Answer the new active
        alarm, if any."""
        time = self.alarmIndex.pop()
        self.alarms.pop(time)
        return self.activeAlarm


    from time import time as getCurrentTime
    from Alarm import Alarm
