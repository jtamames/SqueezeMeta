# AMOS message reader/writer
# Contributed by Paul Harrison

import re

class Error(Exception): pass

class Message:
    """ AMOS Message object
    
        Fields:
	    type         - message type code
	    fields       - dictionary of fields
	    sub_messages - list of sub-messages
	    
	repr(message) converts the message back to AMOS format
    """

    def __init__(self, type):
        self.type = type
	self.fields = { }
	self.sub_messages = [ ]
    
    def __repr__(self):
        result = '{' + self.type + '\n'
	for key in self.fields:
	    result += key+':'
	    if '\n' in self.fields[key]:
	        result += '\n'+self.fields[key]
		if not self.fields[key].endswith('\n'):
		    result += '\n'
		result += '.\n'
	    else:
	        result += self.fields[key]+'\n'	       
	result += ''.join([ repr(sub_message) for sub_message in self.sub_messages ]) + \
	          '}\n'
	return result

    def get_multiline(self, name):
        """ Strip newline characters from a multi-line field. """
        return ''.join(self.fields[name].split('\n'))
    
    def set_multiline(self, name, data):
        """ Set a multi-line field, inserting newlines every 60 characters. """
        self.fields[name] = '\n'.join([
	    data[i:i+60] for i in xrange(0,len(data),60)
	]) + '\n'	       


_START = re.compile(r'^{([A-Z][A-Z][A-Z])\n$')
_MULTILINE_FIELD = re.compile(r'^([a-z][a-z][a-z]):\n$')
_FIELD = re.compile(r'^([a-z][a-z][a-z]):(.*)\n$')

def read_record(file, first_line=None):
    """ Read a record from a file of AMOS messages 
    
        On success returns a Message object	
	On end of file raises EOFError	
	On syntax error raises amos.Error
    """

    if first_line is None:
        first_line = file.readline()
    
    if not first_line:
        raise EOFError()
    
    match = _START.match(first_line)
    if not match:
        raise Error('Bad start of message', first_line)
    
    message = Message(match.group(1))
    
    while True:
        line = file.readline()
	
	match = _MULTILINE_FIELD.match(line)
	if match:
	    name = match.group(1)
	    message.fields[name] = ''
	    while True:
		line = file.readline()
	        if line == '.\n': break
                message.fields[name] += line
	    continue

	match = _FIELD.match(line)
	if match:
	    message.fields[match.group(1)] = match.group(2)
	    continue

	if line == '}\n':
	    break
    
        if line.startswith('{'):
            message.sub_messages.append( read_record(file, line) )
	    continue
	
	raise Error('Bad line',line)

    return message

def iter_records(file):
    """ Iterate over all the records in a file """

    while True:
        try:
            yield read_record(file)
        except EOFError:
	    break


if __name__ == '__main__':  
    #Example: pass a message stream from stdin to stdout
    #         print count of each message type to stderr
    
    import sys
    
    counts = { }
    for record in iter_records(sys.stdin):
        counts[record.type] = counts.get(record.type,0) + 1
        sys.stdout.write(repr(record))
    
    types = counts.keys()
    types.sort()
    for type in types:
        print >> sys.stderr, '%s x %d' % (type, counts[type])
