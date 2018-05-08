# getopt_dict.py
#

"""
Generic helper to wrap the result of getopt(3) into an options dictionary
which is more convenient. See getopt_dict() for help.
"""

import getopt
import sys
import itertools as it

def make_usage(usage_string, retval=1):
    """Given a string, return a usage function which will output the string
    and then terminate the program with the given return code. The returned
    function accepts an optional error string; if given at call time, will
    write the error string to stderr before terminating."""
    def usage_func(uerrstr=''):
        import sys
        sys.stderr.write(usage_string)
        if uerrstr:
            sys.stderr.write('\nError: ' + uerrstr + '\n')
        sys.exit(retval)
    return usage_func

class option(object):
    """Wrapper object for a short or long option as accepted b getopt_dict().
    This is like a POSIX short/long option but supports a special syntax for
    specifying whether the option may repeat (see the getopt_dict() docstring
    for details).

    We internally support optional arguments, but since they are unimplemented
    in Python's getopt module they are currently treated the same as required
    arguments.
    """
    __slots__ = ('_name', '_args', '_repeat', '_value',
            '_start_value', '_short')
    ARG_NONE  = 0
    ARG_YES   = 1
    ARG_MAYBE = 2
    def __init__(self, optstr, short=False):
        optstr = optstr.lstrip('-')
        self._short = short
        self._repeat = False
        # Repeat spec
        if optstr.endswith('*'):
            self._repeat = True
            optstr = optstr[:-1]
        # Core option name
        self._name = optstr.rstrip(':=')
        # Required/optional/no args is determined by trailing ':' or '='
        self._args = len(optstr) - len(self._name) # ARG_*: 0, 1, or 2
        # Base value depends on repeat spec and arg type
        if self._repeat:
            if self._args >= option.ARG_YES:
                self._value = []
            else:
                self._value = 0
        else:
            self._value = False
        # Remember starting value so we know quickly when we have been asserted.
        self._start_value = self._value

    def set(self, arg=None):
        """Indicate that the option has been set on the parsed getopt line.
        This will update the value in various ways depending on whether the
        option can repeat and whether it expects (or requires) an argument.
        An argument of None to this method means no argument is given, but the
        flag is asserted.
        If it is an error to give the specified argument value to this option,
        getopt.GetoptError will be raised."""
        if arg is not None and self._args == option.ARG_NONE:
            self.error("option '%s' does not accept arguments")
        elif arg is None and self._args == option.ARG_YES:
            self.error("option '%s' requires an argument")
        elif arg is None and self._args == option.ARG_MAYBE:
            arg = True
        if self._repeat:
            if self._args >= option.ARG_YES:
                self._value.append(arg)
            else:
                self._value += 1
        else:
            if self.present:
                self.error("option '%s' given multiple times")
            if self._args >= option.ARG_YES:
                self._value = arg
            else:
                self._value = True

    def combine(self, other):
        """Combine this option with another as if they were the same option
        and return the union of the two. It is an error to combine two present
        options which cannot repeat."""
        if other is not None and isinstance(other, option):
            if self._repeat:
                if self._args >= option.ARG_YES:
                    self._value.extend(other._value)
                else:
                    self._value += other._value
            else:
                if self.present and other.present:
                    self.error("option '%%s' given multiple times (as '%s%s')"
                            % (other.prefix(), other.name))
                if other.present:
                    self._value = other._value
        # Take the long option's name if they are different.
        if self.short and other.long:
            self._name = other._name
        return self

    def error(self, fmt_msg):
        """Shortcut for raising getopt.GetoptError."""
        namestr = self.prefix() + self._name
        raise getopt.error(fmt_msg % (namestr,))

    @property
    def value(self):
        """The value of this option, as described in getopt_dict()."""
        return self._value
    @property
    def name(self):
        """The core name of this option, without any formatting characters."""
        return self._name
    @property
    def repeat(self):
        """Whether this option may be specified multiple times (bool)."""
        return self._repeat
    @property
    def present(self):
        """Whether this option has been asserted (is present)."""
        return self._value != self._start_value
    @property
    def short(self):
        """True iff this is a short option (e.g. -x)."""
        return self._short
    @property
    def long(self):
        """True iff this is a long option (e.g. --xx)."""
        return not self._short
    def prefix(self):
        """Return a dash-string based on option type ('-' or '--')."""
        return ('-' if self._short else '--')
    def encode(self):
        """Return an option string suitable for passing to getopt()."""
        optchar = ':' if self._short else '='
        # XXX python does not support '::' or '==', so do not use double chars
        #return self.name + (optchar*self._args)
        return self.name + (optchar if self._args else '')
    def __str__(self):
        return self.encode() + ('*' if self._repeat else '')
    def __repr__(self):
        return '{' + str(self) \
                + (('>'+repr(self.value)) if self.present else '')\
                + '}'

def getopt_dict(args, sopts, lopts=None, stol=dict()):
    """Given a program arguments sequence and short and long option sequences,
    map the given arguments to a dictionary keyed by option names with no
    leading or trailing characters (such as dashes), and with values according
    to whether the option is supplied or not (and how many times - see below).
    Return a tuple (options_dictionary, arguments) similar to getopt().

    If stol is given, it is a dictionary which maps short to long option names,
    in which case the returned dictionary will contain _only_ the long option
    name for those options named in stol. For convenience, the keys/values of
    this dictionary may contain the same formatting as sopts and lopts
    themselves. For example, all short options could be translated into their
    long equivalents for a parallel list pair; stol=dict(zip(sopts, lopts)).

    The format of the options is as with the 'getopt' module (as with POSIX
    getopt) with one extension: if an option ends in *, it may be
    repeated multiple times (otherwise repetition is an error).
    For such an option, the dictionary value will be
    the number of times it is repeated 9for options which do not expect
    arguments) or a list of arguments (for options that accept arguments).
    If an option does not accept repeats, the dictionary value for that option
    will be True/False (for options which do not expect arguments) or the given
    argument/None (for options that accept arguments) when present/not present.


    Note that while this function has internal support for it, the Python
    'getopt' module (at least as late as 2.7.X) does not support optional
    arguments ('::' or '==' suffix in the option string). That means until
    support for optional arguments is implemented, it is not an error to
    provide the '::' or '==' suffix; however, the option will be interpreted as
    always requiring an argument (equivalently to ':' or '=').
    """
    # Create dictionary of special option objects which will handle flag logic.
    options = [option(sopt, True) for sopt in sopts] \
            + [option(lopt, False) for lopt in (lopts or [])]
    options = dict((o.name, o) for o in options)
    # Encode options as string suitable for getopt() (e.g. remove '*' and '-').
    sstrs = ''
    lstrs = []
    for opt in options.values():
        if opt.short:
            sstrs += opt.encode()
        else:
            lstrs.append(opt.encode())
    # Invoke getopt() and update the option dictionary.
    opts, args = getopt.getopt(args, sstrs, lstrs)
    for optopt, optarg in opts:
        # Nb. getopt returns strings as -X, --XX, or --XX= (but never -X:)
        stripped = optopt.lstrip('-')
        opt = options[stripped]
        opt.set(optarg or None)
    # Map the result through the short-to-long dictionary.
    stol = dict(stol)
    for sname, lname in stol.items():
        # Strip names before doing the key lookup.
        sname = option(sname).name
        lname = option(lname).name
        if sname in options:
            options[lname] = options[sname].combine(options.get(lname))
            del options[sname]
    # Normalize the return dictionary by inserting final option values
    # (but only for options which are present).
    optdict = dict((stol.get(opt.name, opt.name), opt.value)
                    for opt in options.values() if opt.present)
    return (optdict, args)

