__all__ = ['atomic_number', 'atomic_mass']

def atomic_number(s: str):
    d = {'h': 1, 'he': 2, 'li': 3, 'be': 4, 'b': 5, 'c': 6, 'n': 7, 'o': 8,
         'f': 9, 'ne': 10, 'na': 11, 'mg': 12, 'al': 13, 'si': 14, 'p': 15,
         's': 16, 'cl': 17, 'ar': 18, 'k': 19, 'ca': 20, 'sc': 21, 'ti': 22,
         'v': 23, 'cr': 24, 'mn': 25, 'fe': 26, 'co': 27, 'ni': 28, 'cu': 29,
         'zn': 30, 'ga': 31, 'ge': 32, 'as': 33, 'se': 34, 'br': 35, 'kr': 36,
         'rb': 37, 'sr': 38, 'y': 39, 'zr': 40, 'nb': 41, 'mo': 42, 'tc': 43,
         'ru': 44, 'rh': 45, 'pd': 46, 'ag': 47, 'cd': 48, 'in': 49, 'sn': 50,
         'sb': 51, 'te': 52, 'i': 53, 'xe': 54, 'cs': 55, 'ba': 56, 'la': 57,
         'ce': 58, 'pr': 59, 'nd': 60, 'pm': 61, 'sm': 62, 'eu': 63, 'gd': 64,
         'tb': 65, 'dy': 66, 'ho': 67, 'er': 68, 'tm': 69, 'yb': 70, 'lu': 71,
         'hf': 72, 'ta': 73, 'w': 74, 're': 75, 'os': 76, 'ir': 77, 'pt': 78,
         'au': 79, 'hg': 80, 'tl': 81, 'pb': 82, 'bi': 83, 'po': 84, 'at': 85,
         'rn': 86, 'fr': 87, 'ra': 88, 'ac': 89, 'th': 90, 'pa': 91, 'u': 92,
         'np': 93, 'pu': 94, 'am': 95, 'cm': 96, 'bk': 97, 'cf': 98, 'es': 99,
         'fm': 100, 'md': 101, 'no': 102, 'lr': 103, 'rf': 104, 'db': 105}
    return d[s]


def atomic_mass(s: str):
    d = {'h': 1.008, 'he': 4.003, 'li': 6.941, 'be': 9.012, 'b': 10.811, 'c': 12.011, 'n': 14.007, 'o': 15.999,
         'f': 18.998, 'ne': 20.180, 'na': 22.990, 'mg': 24.305, 'al': 26.982, 'si': 28.086, 'p': 30.974,
         's': 32.065, 'cl': 35.453, 'ar': 39.948, 'k': 39.098, 'ca': 40.078, 'sc': 44.956, 'ti': 47.867,
         'v': 50.942, 'cr': 51.996, 'mn': 54.938, 'fe': 55.845, 'co': 58.933, 'ni': 58.693, 'cu': 63.546,
         'zn': 65.390, 'ga': 69.723, 'ge': 72.640, 'as': 74.922, 'se': 78.960, 'br': 79.904, 'kr': 83.800,
         'rb': 85.468, 'sr': 87.620, 'y': 88.906, 'zr': 91.224, 'nb': 92.906, 'mo': 95.940, 'tc': 98.000,
         'ru': 101.070, 'rh': 102.906, 'pd': 106.420, 'ag': 107.868, 'cd': 112.411, 'in': 114.818, 'sn': 118.710,
         'sb': 121.760, 'te': 127.600, 'i': 126.905, 'xe': 131.293, 'cs': 132.906, 'ba': 137.327, 'la': 138.906,
         'ce': 140.116, 'pr': 140.908, 'nd': 144.240, 'pm': 145.000, 'sm': 150.360, 'eu': 151.964, 'gd': 157.250,
         'tb': 158.925, 'dy': 162.500, 'ho': 164.930, 'er': 167.259, 'tm': 168.934, 'yb': 173.040, 'lu': 174.967,
         'hf': 178.490, 'ta': 180.948, 'w': 183.840, 're': 186.207, 'os': 190.230, 'ir': 192.217, 'pt': 195.078,
         'au': 196.967, 'hg': 200.590, 'tl': 204.383, 'pb': 207.200, 'bi': 208.980, 'po': 209.000, 'at': 210.000,
         'rn': 222.000, 'fr': 223.000, 'ra': 226.000, 'ac': 227.000, 'th': 232.038, 'pa': 231.036, 'u': 238.029,
         'np': 237.000, 'pu': 244.000, 'am': 243.000, 'cm': 247.000, 'bk': 247.000, 'cf': 251.000, 'es': 252.000,
         'fm': 257.000, 'md': 258.000, 'no': 259.000, 'lr': 262.000, 'rf': 261.000, 'db': 262.000}
    return d[s]

