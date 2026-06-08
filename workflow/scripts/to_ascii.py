#!/usr/bin/env python3
"""
Transliterate non-ASCII characters to ASCII equivalents (stdin -> stdout).

Handles Greek letters explicitly; falls back to NFKD decomposition for
accented Latin characters (e.g. e-acute -> e). Remaining unmappable
characters are replaced with an underscore.

Drop-in replacement for `iconv -f UTF-8 -t ASCII//TRANSLIT`, which is
locale-dependent and silently produces '?' for Greek in LANG=C environments.
"""
import sys
import unicodedata

_GREEK = dict(zip(
    'αβγδεζηθικλμνξοπρστυφχψωΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ',
    ('alpha beta gamma delta epsilon zeta eta theta iota kappa lambda mu nu '
     'xi omicron pi rho sigma tau upsilon phi chi psi omega '
     'Alpha Beta Gamma Delta Epsilon Zeta Eta Theta Iota Kappa Lambda Mu Nu '
     'Xi Omicron Pi Rho Sigma Tau Upsilon Phi Chi Psi Omega').split()
))

def to_ascii(text):
    out = []
    for c in text:
        if ord(c) < 128:
            out.append(c)
        elif c in _GREEK:
            out.append(_GREEK[c])
        else:
            norm = unicodedata.normalize('NFKD', c).encode('ascii', 'ignore').decode('ascii')
            out.append(norm if norm else '_')
    return ''.join(out)

for line in sys.stdin:
    sys.stdout.write(to_ascii(line))
