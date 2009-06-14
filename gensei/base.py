import os, sys
import numpy as np
import copy
import types

from getch import getch

def import_file(filename):
    path = os.path.dirname(filename)
    if not path in sys.path:
        sys.path.append(path)
    name = os.path.splitext(os.path.basename(filename))[0]
    mod = __import__(name)
    return mod

def assert_(condition):
    if not condition:
        raise ValueError('assertion failed!')

def pause( msg = None ):
    """
    Prints the line number and waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging.
    """
    f = sys._getframe(1)
    ff = f.f_code
    if (msg):
        print '%s, %d: %s(), %d: %s' % (ff.co_filename, ff.co_firstlineno,
                                        ff.co_name, f.f_lineno, msg)
    else:
        print '%s, %d: %s(), %d' % (ff.co_filename, ff.co_firstlineno,
                                    ff.co_name, f.f_lineno)
    spause()

def spause( msg = None ):
    """
    Waits for a keypress.

    If you press:
    "q" ............. it will call sys.exit()
    any other key ... it will continue execution of the program

    This is useful for debugging. This function is called from pause().
    """
    if (msg):
        print msg
    sys.stdout.flush()
    ch = getch()
    if ch == 'q':
        sys.exit()

class Object(object):
    traits = {}

    @staticmethod
    def objects_from_dict(*args, **kwargs):

        kws = copy.copy(kwargs)
        
        level = kws.pop('level', 0)
        flag = kws.pop('flag')

        # For level 0 only...
        constructor = kws.pop('constructor', Object)

        out = []
        for arg in args:
            if type(arg) == dict:
                if flag[level]:
                    arg.update(kws)
                    aux = constructor(**arg)
                    iterator = aux.__dict__
                else:
                    aux = iterator = arg

                for key, val in iterator.iteritems():
                    if (type(val) == dict) and (key != 'traits'):
                        try:
                            flag[level+1]
                        except:
                            flag = flag + (0,)
                        val2 = Object.objects_from_dict(val, level=level+1,
                                                        flag=flag)

                        if flag[level]:
                            aux.__dict__[key] = val2
                        else:
                            aux[key] = val2

                out.append(aux)
            else:
                out.append(arg)

        if len(out) == 1:
            out = out[0]

        return out

    def __init__(self, **kwargs):
        if kwargs:
            self.__dict__.update(kwargs)
            self.traits = copy.copy(self.__class__.traits)
            self.set_default_traits(kwargs.iterkeys())

    def set_default_traits(self, keys):
        self.traits.update({}.fromkeys(keys))

    def init_trait(self, name, value, format=None):
        setattr(self, name, value)
        self.traits[name] = format
            
    def __str__(self):
        msg = ['%s' % object.__str__(self)]

        keys = self.traits.keys()
        order = np.argsort(keys)
        for ii in order:
            key = keys[ii]
            val = self.traits[key]

            if isinstance(val, tuple):
                tr = val[1]
                try:
                    attr = tr(getattr(self, key))
                except:
                    attr = '<not set>'
                    val = '%s'
                else:
                    val = val[0]
            elif isinstance(val, types.FunctionType):
                attr = val(self)
                val = None
            else:
                try:
                    attr = getattr(self, key)
                except:
                    attr = '<not set>'
                
            if issubclass(attr.__class__, Object):
                sattr = repr(attr)
                attr = '%s: %s' % (key, sattr)
            else:
                if val is None:
                    attr = '%s: %s' % (key, attr)
                else:
                    attr = ('%s: ' % key) + (val % attr)

            msg.append(attr)

        return '\n'.join(msg)

    def __repr__(self):
        if hasattr(self, 'name'):
            return "%s:%s" % (self.__class__.__name__, self.name)
        else:
            return object.__repr__(self)

    def copy(self, deep=False):
        if deep:
            return copy.deepcopy(self)
        else:
            return copy(self)

class Config(Object):

    @staticmethod
    def from_file(filename, required, optional):
        conf_mod = import_file(filename)

        if 'define' in conf_mod.__dict__:
            define_dict = conf_mod.__dict__['define']()
        else:
            define_dict = conf_mod.__dict__

        return Config.from_conf(define_dict, required, optional)

    @staticmethod
    def from_conf(conf, required, optional):
        valid = {}
        for kw in required:
            try:
                val = conf[kw]
            except KeyError:
                raise ValueError('missing keyword "%s" in "%s"!'
                                 % (kw, filename))
            valid[kw] = val

        for kw in optional:
            valid[kw] = conf.get(kw, None)

        return Config(**valid)

    def override(self, options, can_override):
        for group_key, group in self.__dict__.iteritems():
            for key, val in group.iteritems():
                if key in can_override:
                    group[key] = getattr(options, key)

class Output(Object):
    """Factory class providing output (print) functions.

    Example:

    >>> output = Output( 'sfepy:' )
    >>> output( 1, 2, 3, 'hello' )
    >>> output.prefix = 'my_cool_app:'
    >>> output( 1, 2, 3, 'hello' )
    """
    traits = {
        'prefix' : None,
        'output_function' : None,
        'level' : None,
    }

    def __init__(self, prefix, filename=None, combined=False, **kwargs):
        Object.__init__(self, **kwargs)

        self.prefix = prefix

        self.set_output(filename, combined)
        
    def __call__(self, *argc, **argv):
        self.output_function(*argc, **argv)

    def set_output(self, filename=None, combined=False, append=False):
        """Set the output function - all SfePy printing is accomplished by
        it. If filename is None, output is to screen only, otherwise it is to
        the specified file, moreover, if combined is True, both the ways are
        used.

        Arguments:
                filename - print into this file
                combined - print both on screen and into a file
                append - append to an existing file instead of overwriting it
        """
        self.level = 0
        def output_screen( *argc, **argv ):
            format = '%s' + ' %s' * (len( argc ) - 1)
            msg =  format % argc

            if msg.startswith( '...' ):
                self.level -= 1

            print self._prefix + ('  ' * self.level) + msg

            if msg.endswith( '...' ):
                self.level += 1

        def output_file( *argc, **argv ):
            format = '%s' + ' %s' * (len( argc ) - 1)
            msg =  format % argc

            if msg.startswith( '...' ):
                self.level -= 1

            fd = open( filename, 'a' )
            print >>fd, self._prefix + ('  ' * self.level) + msg
            fd.close()

            if msg.endswith( '...' ):
                self.level += 1

        def output_combined( *argc, **argv ):
            output_screen( *argc, **argv )
            output_file( *argc, **argv )
    
        if filename is None:
            self.output_function = output_screen

        else:
            if not append:
                fd = open( filename, 'w' )
                fd.close()

            if combined:
                self.output_function = output_combined
            else:
                self.output_function = output_file

    def get_output_function(self):
        return self.output_function

    def set_output_prefix( self, prefix ):
        assert_( isinstance( prefix, str ) )
        if len( prefix ) > 0:
            prefix += ' '
        self._prefix = prefix
        
    def get_output_prefix( self ):
        return self._prefix[:-1]
    prefix = property( get_output_prefix, set_output_prefix )
    
output = Output('gensei:')
