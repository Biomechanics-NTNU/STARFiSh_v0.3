
__author__ = "Vinzenz Gregor Eck"
__version__ = "0.3"

def systemCheck():
    '''
    Check for necessary moduls and 3thrd party libarys needed for STARFiSh and vnc:
        matplotlib
        numpy
        scipy
        lxml2
        pydot
        gtk
        pyglet
        mencoder
        OpenGL
    '''
    matplotlibInstalled = False
    installed = " ... installed!"
    print ""
    print "Check for necessary moduls and 3thrd party libarys needed for vascular1DFlow and vnc"
    print ""
    
    print "check for matplotlib ",
    try:
        import matplotlib.pyplot 
        print installed
        matplotlibInstalled = True
    except:
        print " IMPORT ERROR; no version of matplotlib.pyplot found"
    
    print "check for mencoder   ",
    try:
        if matplotlibInstalled:
            import matplotlib.animation as animation
            
            if animation.MencoderWriter.isAvailable():
                print installed
            else:
                print " IMPORT ERROR; no version of mencoder found"     
        else:
            print " IMPORT ERROR; no version of matplotlib.pyplot found thus cannot check for mencoder"     
    except:    
        print " IMPORT ERROR; no version of mencoder found"    
        
    print "check for numpy      ",
    try:
        import numpy
        print installed
        
    except:
        print " IMPORT ERROR; no version of numpy found"
    
    print "check for scipy      ",
    try:
        import scipy
        print installed
    except:
        print " IMPORT ERROR; no version of scipy found"
        
    print "check for lxml2      ",
    try:
        import lxml
        print installed
    except:
        print " IMPORT ERROR; no version of lxml2 found"
    
    print "check for pydot      ",
    try:
        import pydot
        print installed
    except:
        print " IMPORT ERROR; no version of pydot found"
    
    print "check for gtk        ",
    try:
        import gtk
        print installed
    except:
        print " IMPORT ERROR; no version of gtk found" 
    
    print "check for pyglet     ",
    try:
        import pyglet
        print installed
    except:
        print " IMPORT ERROR; no version of pyglet found" 
        
    print "check for pyOpenGL   ",
    try:
        import OpenGL
        print installed
    except:
        print " IMPORT ERROR; no version of pyOpenGL found" 
        
    print "check for pydot      ",
    try:
        import pydot
        print installed
    except:
        print "IMPORT ERROR; no version of found (http://pypi.python.org/pypi/pydot/1.0.2)" 


if __name__ == '__main__':
    systemCheck()