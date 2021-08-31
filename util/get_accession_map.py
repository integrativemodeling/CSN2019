# Create the accession map used by test/mock/sitecustomize.py

import ihm.reference

codes = ['Q13098', 'P61201', 'Q9UNS2', 'Q9BT78', 'Q92905', 'Q7L5N1',
         'Q9H9Q2', 'Q99627', 'Q8WXC6']

def pp(s):
    indent = 8
    width = 66
    def get_lines(s):
        for i in range(0, len(s), width):
            yield ' ' * indent + "'" + s[i:i+width] + "'"
    return '\n'.join(l for l in get_lines(s))

for code in codes:
    u = ihm.reference.UniProtSequence.from_accession(code)
    print("    '%s': {'db_code':'%s', 'sequence':\n%s},"
          % (code, u.db_code, pp(u.sequence)))
