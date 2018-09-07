## {{{ http://code.activestate.com/recipes/577871/ (r1) + modification by Alessio Porcelli (ale_led@yahoo.it)
import sys
import time

class ProgressBar(object):
    """ProgressBar class holds the options of the progress bar.
    The options are:
        start   State from which start the progress. For example, if start is 
                5 and the end is 10, the progress of this state is 50%
        end     State in which the progress has terminated.
        width   --
        fill    String to use for "filled" used to represent the progress
        blank   String to use for "filled" used to represent remaining space.
        format  Format
        incremental
    """
    def __init__(self, start=0, end=10, width=12, fill='=', blank='.', format='[{fill}>{blank}] {progress}%', incremental=True):
        super(ProgressBar, self).__init__()
        self.start = start
        self.end = end
        self.width = width
        self.fill = fill
        self.blank = blank
        self.format = format
        self.incremental = incremental
        self.step = 100 / float(width) #fix
        self.reset()

    def __add__(self, increment):
        increment = self._get_progress(increment)
        #if 100 > self.progress + increment:
        self.progress += increment
        #else:
        #    self.progress = 100
        return self

    def __str__(self):
        progressed = int(round(self.progress / self.step,0)) #fix
        fill = progressed * self.fill
        blank = (self.width - progressed) * self.blank
        return self.format.format(fill=fill, blank=blank, progress=int(round(self.progress,0)))

    __repr__ = __str__

    def _get_progress(self, increment):
        return float(increment * 100) / self.end

    def reset(self):
        """Resets the current progress to the start point"""
        self.progress = self._get_progress(self.start)
        return self


class AnimatedProgressBar(ProgressBar):
    """Extends ProgressBar to allow you to use it straighforward on a script.
    Accepts an extra keyword argument named `stdout` (by default use sys.stdout)
    and may be any file-object to which send the progress status.
    """
    def __init__(self, *args, **kwargs):
        self.stdout = kwargs.pop('stdout', sys.stdout)
        super(AnimatedProgressBar, self).__init__(*args, **kwargs)

    def show_progress(self, info=None):
        if hasattr(self.stdout, 'isatty') and self.stdout.isatty():
            # cariage return + kill the line! ...ready for a new one!
            #\33[ is ESC; K is kill: K or 0K form cursor to end, 1K from cursor to begin, 2K all the line. Cursor position does not change.
            self.stdout.write('\r\33[2K')
        else:
            self.stdout.write('\e[2K\n') # How to kill the line in non ANSI, as windows? Is it correct? UNTESTED!
        self.stdout.write(str(self)+("" if info is None else "  %s"%info))
        self.stdout.flush()

    def progressing(self, increment, info=None):
        self + increment
        self.show_progress(info)

    def finish(self, info=None):
        self.show_progress(("" if info is None else info)+'\n')


if __name__ == '__main__':
    p = AnimatedProgressBar(end=100, width=80)

    while True:
        p + 5
        p.show_progress()
        time.sleep(0.1)
        if p.progress == 100:
            break
    print #new line
    print 'OR'
    while True:
        p.progressing(5,'just a test...') # p+5 and show_progress in one shot. Moreover an info is printed
        time.sleep(0.1)
        if p.progress == 100:
            break
    p.finish("it's over.") # end (and an info is printed)
## end of http://code.activestate.com/recipes/577871/ + modification by Alessio Porcelli (ale_led@yahoo.it) }}}
