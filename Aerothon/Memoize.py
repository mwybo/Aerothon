"""
A decorator for memoizing functions

http://pko.ch/2008/08/22/memoization-in-python-easier-than-what-it-should-be/

"""
#import functools
#import cPickle
#def memoize(fctn):
#        memory = {}
#        @functools.wraps(fctn)
#        def memo(*args,**kwargs):
#                haxh = cPickle.dumps((args, sorted(kwargs.iteritems())))
#
#                if haxh not in memory:
#                        memory[haxh] = fctn(*args,**kwargs)
#
#                return memory[haxh]
#        if memo.__doc__:
#            memo.__doc__ = "\n".join([memo.__doc__,"This function is memoized."])
#        return memo


class memoize:
    """
    http://snippets.dzone.com/posts/show/4840
    """
    
    def __init__ (self, f):
        self.f = f
        self.mem = {}
    def __call__ (self, *args, **kwargs):
        if (args, str(kwargs)) in self.mem:
            print "Used memoization!"
            return self.mem[args, str(kwargs)]
        else:
            tmp = self.f(*args, **kwargs)
            self.mem[args, str(kwargs)] = tmp
            return tmp
