


class AlignedPrinter:
  def __init__(self):
    self._lefts = []
    self._rights = []

  def add(self, left, right):
    self._lefts.append(left)
    self._rights.append(right)

  def sort_right_float(self):
    floating_array = []
    for elem in self._rights:
      floating_array.append(float(elem.split()[0]))
    self._lefts = [x for _,x in sorted(zip(floating_array, self._lefts))]
    self._rights = [x for _,x in sorted(zip(floating_array, self._rights))]

  def display(self):
    max_chars = 0
    for l in self._lefts:
      max_chars = max(max_chars, len(l))
    for i in range(0, len(self._lefts)):
      l = self._lefts[i]
      r = self._rights[i]
      to_print = l + " " * (max_chars - len(l) + 1) + r
      print(to_print)

