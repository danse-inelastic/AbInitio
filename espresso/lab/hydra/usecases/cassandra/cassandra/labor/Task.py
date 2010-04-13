

class Task(object):


    def __init__(self, id, input):
        self.id = id
        self.input = input
        return


    def __repr__(self):
        return "<task %r>" % self.id