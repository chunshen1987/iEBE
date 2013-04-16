#! /usr/bin/env python
# This script translates a list of arguments into one value specified by a rule file. Usage:
# translate ruleFilename key1 key2 ...
# It prints out the value corresponds to [key1, key2, ...] from a dictionary read from ruleFilename.
# To see how the dictionary is generated, see readRules function.

from sys import exit, argv

def processOneLine(aLine, level_indicator="+", key_separator=":", commentSymbol="#"):
  """
    Return [level, keys_list, value] list from string aLine.
    level is indicated by how many successive level_indicators are there to the left, key and value are separated by key_separator.
  """
  # take care of comments:
  if commentSymbol in aLine:
    aLine = aLine[:aLine.index(commentSymbol)].strip();

  # if it's an empty line:
  aLine = aLine.strip()
  if aLine=="": return []

  # check if syntax is correct:
  if key_separator not in aLine:
    print("translate.processOneLine error: key-value separator "+key_separator+" not included in the line \n"+aLine)
    exit(-1)

  # get level
  level = 0
  for i in range(len(aLine)):
    if aLine[i]==level_indicator:
      level = level + 1
    else:
      aLine = aLine[i:]
      break

  # separate key and value
  components = aLine.split(key_separator);
  keys_list = [x.strip() for x in components[:-1]];
  value = components[-1].strip();

  # finally...
  return [level, keys_list, value]


def readRules(buffer,level_indicator="+", key_separator=":", commentSymbol="#"):
  """
    Process the text buffer to get the rule used for translations, line by line. Each line will be transferred into one entry in a rule dictionary. The dictionary will then be returned. The rule dictionary is generated using all the list of all strings between key_separators except the last one as the key, and the last one as value.
    For example,
    a : b: 1 # comments
    will be translates into entry ["a","b"]:"1"
    To ease the pain for repeated common keys, a level_indicator can be used to indicate how may shared keys the current line inherits from previous lines. The number of level_indicator means the number of keys the current line should inherit (starts from left) from previous lines.
    For example, if the text buffer looks like:
    z : 1
    + a : 2
    ++ b : 3
    + d: 4
    The rule dictionary will contain:
    ("z") : "1"
    ("z", "a") : 2
    ("z", "a", "b") : 3
    ("z", "d") : 4
    Note that the following
    z : 1
    ++ a : 2
    will raise an error.
  """
  D = {}
  accumulated_keys = [];
  for aLine in buffer:
    tmp_result = processOneLine(aLine)
    if not tmp_result: continue
    level, keys, value = tmp_result
    if level>len(accumulated_keys):
      print("translates.readRules error: two many "+level_indicator+" signs in the line\n"+aLine)
      exit(-1)
    else:
      accumulated_keys = accumulated_keys[:level]
      accumulated_keys.extend(keys)
      D[tuple(accumulated_keys)] = value
  return D


def translate(ruleFilename, keys_list):
  """
    Translate keys_list into the correponding value given in the dictionary generated from ruleFilename using readRules function.
  """
  D = readRules(ruleFilename)
  result = ""
  for ii in range(len(keys_list)): result+=" "+(D[tuple(keys_list[:ii+1])])
  return result

if __name__ == '__main__':
  if len(argv)<3:
    print("Usage: translate ruleFilename key1 key2 ...")
    exit(-1)
  else:
    print(translate(file(argv[1]).readlines(),argv[2:]))
