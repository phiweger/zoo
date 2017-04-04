import re


# valid
a = 'A/mallard/Interior Alaska/6MP0155/2006(H3N8)'
b = 'A/blue-winged Teal/Minnesota/AI09-2977/2009(H4N8)'
c = 'A/Peru/PER175/2012(H3N2)'

# not valid

d = 'A/mallard/Interior Alaska/6MP0155/2006'
e = 'A/Peru/PER175/4912(H3N2)'


# stackoverflow, 4374185
p = re.compile('A/[.+?/]?.+?/.+?/(19|20)\d{2}\(H\dN\d\)')


# match
print(re.match(p, a).group(0))
print(re.match(p, b).group(0))
print(re.match(p, c).group(0))

# no match, AttributeError: 'NoneType' object has no attribute 'group'
print(re.match(p, d).group(0))
print(re.match(p, e).group(0))
