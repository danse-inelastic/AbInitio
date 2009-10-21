#!/usr/bin/env python

# A very simple example of python curses
# Taken from: heather.cs.ucdavis.edu/~matloff/Python/PyCurses.pdf
# Usage:    python testcrs.py <number>
# Example:  python testcrs.py 5


import curses, sys, traceback
class gb:
  boxrows = int(sys.argv[1]) # number of rows in the box
  boxcols = boxrows # number of columns in the box
  scrn = None # will point to window object
  row = None # current row position
  col = None # current column position
  
def draw(chr):
  if gb.row == gb.boxrows-1:
     gb.scrn.addch(gb.row,gb.col,chr,curses.color_pair(1))
  else:
     gb.scrn.addch(gb.row,gb.col,chr)
  gb.scrn.refresh()
  gb.row += 1
  if gb.row == gb.boxrows:
     gb.row = 0
     gb.col += 1
     if gb.col == gb.boxcols: gb.col = 0

def restorescreen():
  curses.nocbreak()
  curses.echo()
  curses.endwin()

def main():
  gb.scrn = curses.initscr()
  curses.noecho()
  curses.cbreak()
  curses.start_color()
  curses.init_pair(1,curses.COLOR_RED,curses.COLOR_WHITE)
  gb.scrn.clear()
  gb.row = 0
  gb.col = 0
  gb.scrn.refresh()
  while True:
     c = gb.scrn.getch()
     c = chr(c)
     if c == 'q': break
     draw(c)
  restorescreen()
  
if __name__ =='__main__':
  try:
     main()
  except:
     restorescreen()
     traceback.print_exc()


__date__ = "$Oct 21, 2009 6:52:46 AM$"


