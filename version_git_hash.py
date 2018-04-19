import subprocess
import os

this_files_directory = os.path.dirname(os.path.abspath(__file__))
completed_obj = subprocess.run(['git', 'rev-parse', 'HEAD'], stdout=subprocess.PIPE, check=True, cwd=this_files_directory)
git_hash = completed_obj.stdout.decode('utf-8').strip()
print(git_hash)
