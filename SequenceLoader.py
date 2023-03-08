class SequenceLoader:
    def __init__(self):
        pass
    def load( self, FILENAME ):
      raw = open(FILENAME,'r').readlines()
      sequences = {}

      i = 0
      cur_id  = ''
      cur_str = ''
      for liner in raw:
        line = liner.strip()
        if '#' in line:
          if len(cur_id) > 0:
            sequences[i] = (i,cur_id,cur_str)
            cur_id  = ''
            cur_str = ''
            i += 1
          parts = line.split('#')
          cur_id = parts[0]
          cur_str += parts[1]
        else:
          cur_str += line
      if len(cur_str) > 0 or len(cur_id) > 0:
        sequences[i] = (i,cur_id,cur_str)
      return sequences