#!/usr/bin/env python3

import argparse
import subprocess
import re
import json
import os, sys
from pprint import pprint



def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs) 


def getEmoji(string):
  if string == 'Not Found':
    return '&#10060;'
  elif string == 'Skipped':
    return '&#9888;'
  else:
    return '&#9989;'
  
def loadDataFromJson(filename):
  try:
    with open(filename) as json_file:
      return json.load(json_file)
  except Exception as e:
    eprint(f"JSON_ERROR: Unable to load {filename}:\n  {e}")
    return None


def checkString(inputString, regex):
  """
  check if inputString ~/regex/
  """
  pattern = re.compile(regex)
  match = pattern.search(inputString)
  if match:
    return True
  return False
  

def checkBinaryVersion(binary, exitcode=None, text=None, out=None, err=None, version=None, params=[]):
  """
  execute a binary. check if exitcode is correct, stdout containts 'out' and stderr contains 'err'
  """
  command = [binary]
  if len(params) > 0:
    command.extend(params)

  try:
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  except Exception as e:
    if opt.verbose:
      eprint(f"Warning: {binary} was not executed: {e}")
    return [False, None]

  combinedProgramOutput = result.stderr.decode('utf-8') + '\n|\n' + result.stdout.decode('utf-8')
  programVersion=None
  errors = 0

  if text is not None:
    if not checkString(combinedProgramOutput, text):
      if opt.verbose:
        eprint(f"[binary] standard output/error, expected {text}, got {combinedProgramOutput}")
      errors += 1

  if out is not None:
    if not checkString(result.stdout.decode('utf-8'), out):
      if opt.verbose:
        eprint(f"[binary] standard output, expected {out}, got {result.stdout.decode('utf-8')}")
      errors += 1

  if err is not None:
    if not checkString(result.stderr.decode('utf-8'), err):
      if opt.verbose:
        eprint(f"[binary] standard error, expected {err}, got {result.stderr.decode('utf-8')}")
      errors += 1

  if exitcode is not None:
    if result.returncode != exitcode:
      if opt.verbose:
        eprint(f"[binary] exit status, expected {exitcode}, got {result.returncode}")
      errors += 1
  
  if version is not None:
    matches = re.findall(version, combinedProgramOutput)
    if matches:
      programVersion = matches[0]
    else:
      programVersion = 'undefined'

  if errors > 0:
    return [False, None]

  return [True, programVersion]

def checkBinaryFromObject(data):
  exitcode, stdout, stderr, versionRe, anyout = None, None, None, None, None
  parameters = []
  
  if 'binary' in data:
    binary = data['binary']
  else:
    binary = package

  if 'exitcode' in data:
    exitcode = data['exitcode']
      
  if 'stdout' in data:
    stdout = data['stdout']
      
  if 'stderr' in data:
    stderr = data['stderr']

  if 'log' in data:
    anyout = data['log']    
  
  if 'version' in data:
    versionRe = data['version']
      
  if 'params' in data:
    parameters = data['params']
      
  status, ver = checkBinaryVersion(binary, exitcode=exitcode, text=anyout, out=stdout, err=stderr, version=versionRe, params=parameters)
  return [status, ver]

if __name__ == "__main__":
  #databaseFile =  os.path.join(os.path.dirname(os.path.realpath(__file__)), "checker.json")

  file_errors = 0
  package_errors = 0
  file_paths = {}
  package_versions = {}

  opt_parser = argparse.ArgumentParser(description='utility to check the presence of binaries and files')
  opt_parser.add_argument('-j', '--jsonfile',
                        help='Binary to be checked in JSON format',
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "checker.json"))
                         
  opt_parser.add_argument('-f', '--filepath', nargs='+',
                        help='Check the existence of this path')
  opt_parser.add_argument('--html',
                        help='Write an HTML summary')
  opt_parser.add_argument('--verbose',
                        help='Enable verbose output', action='store_true')
  
  opt = opt_parser.parse_args()                     
  data = loadDataFromJson(opt.jsonfile)

  if data is None:
    eprint(f"ERROR: No data loaded from <{opt.jsonfile}>")
    exit()

   
  if 'binaries' in data:
    for package in data['binaries']:
      if opt.verbose:
        eprint('- Checking package', package)
      
      
      status, ver = checkBinaryFromObject(data['binaries'][package])
      if status:
        if opt.verbose:
          eprint(f"\tOK {ver}")
        package_versions[package] = ver
      else:
        if 'optional' in data['binaries'][package] and data['binaries'][package]['optional']:
          if opt.verbose:
            eprint("\tNOT_FOUND_IGNORED")
          package_versions[package] = 'Skipped'
        else:
          if opt.verbose:
            eprint("\tNOT_FOUND_ERROR")
          package_versions[package] = 'Not Found'
          package_errors += 1
  else:
    eprint(f"WARNING: <binaries> section not found in {opt.jsonfile}")

   
  if 'files' in data:
    for fileID in data['files']:
      file = data['files'][fileID]
      if opt.verbose:
        eprint('- Checking file', fileID)

      if os.path.exists(file['path']):
        if opt.verbose:
          eprint("\tOK")
        file_paths[fileID] = "OK"
      elif file['optional']:
        if opt.verbose:
          eprint("\tNOT_FOUND_IGNORED")
        file_paths[fileID] = "Skipped"
      else:
        if opt.verbose:
          eprint("\tNOT_FOUND_ERROR")
        file_paths[fileID] = "Not Found"
        file_errors +=1

  else:
    eprint(f"WARNING: <files> section not found in {opt.jsonfile}")
  
  if opt.filepath:
    for filePath in opt.filepath:
      if opt.verbose:
        eprint('- Checking custom file', filePath)

      if os.path.exists(filePath):
        if opt.verbose:
          eprint("\tOK")
        file_paths[filePath] = "OK"
      else:
        if opt.verbose:
          eprint("\tNOT_FOUND_ERROR")
        file_paths[filePath] = "Not Found"
        file_errors +=1
     

  if len(package_versions) > 0:
    html_string = '<h3>Tested tools</h3>\n<table>'
  for tool in package_versions:
    print(f"Tool:{tool}\t{package_versions[tool]}")
    html_string += f'<tr>\n<td>{getEmoji(package_versions[tool])}</td><td>{tool}</td><td>{package_versions[tool]}</td>\n</tr>\n'
  if len(package_versions) > 0:
    html_string += '</table>'

  if len(file_paths) > 0:
    html_string += '\n<h3>Tested paths</h3>\n<table>\n'
  for file in file_paths:
    print(f"File:{file}\t{file_paths[file]}")
    html_string += f'<tr>\n<td>{getEmoji(file_paths[file])}</td><td>{file}</td><td>{file_paths[file]}</td>\n</tr>\n'
  if len(file_paths) > 0:
    html_string += '</table>\n'

  if opt.html is not None:
    try:
      f = open(opt.html, "w")
      f.write(html_string)
      f.close()
    except Exception as e:
      eprint(f"Error: unable to write HTML to {opt.html}:\n  {e}")

  if package_errors > 0 or file_errors > 0:
    eprint(f"FAILED VALIDATION: {package_errors} packages not found; {file_errors} not found.")