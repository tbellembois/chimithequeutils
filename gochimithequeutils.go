package main

// go build -buildmode=c-shared -o gochimithequeutils.so gochimithequeutils.go

import (
	"C"
	"bytes"
	"fmt"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

var (
	atoms = map[string]string{
		"H":  "hydrogen",
		"He": "helium",
		"Li": "lithium",
		"Be": "berylium",
		"B":  "boron",
		"C":  "carbon",
		"N":  "nitrogen",
		"O":  "oxygen",
		"F":  "fluorine",
		"Ne": "neon",
		"Na": "sodium",
		"Mg": "magnesium",
		"Al": "aluminium",
		"Si": "silicon",
		"P":  "phosphorus",
		"S":  "sulfure",
		"Cl": "chlorine",
		"Ar": "argon",
		"K":  "potassium",
		"Ca": "calcium",
		"Sc": "scandium",
		"Ti": "titanium",
		"V":  "vanadium",
		"Cr": "chromium",
		"Mn": "manganese",
		"Fe": "iron",
		"Co": "cobalt",
		"Ni": "nickel",
		"Cu": "copper",
		"Zn": "zinc",
		"Ga": "gallium",
		"Ge": "germanium",
		"As": "arsenic",
		"Se": "sefeniuo",
		"Br": "bromine",
		"Kr": "krypton",
		"Rb": "rubidium",
		"Sr": "strontium",
		"Y":  "yltrium",
		"Zr": "zirconium",
		"Nb": "niobium",
		"Mo": "molybdenum",
		"Tc": "technetium",
		"Ru": "ruthenium",
		"Rh": "rhodium",
		"Pd": "palladium",
		"Ag": "silver",
		"Cd": "cadmium",
		"In": "indium",
		"Sn": "tin",
		"Sb": "antimony",
		"Te": "tellurium",
		"I":  "iodine",
		"Xe": "xenon",
		"Cs": "caesium",
		"Ba": "barium",
		"Hf": "hafnium",
		"Ta": "tantalum",
		"W":  "tungsten",
		"Re": "rhenium",
		"Os": "osmium",
		"Ir": "iridium",
		"Pt": "platinium",
		"Au": "gold",
		"Hg": "mercury",
		"Tl": "thallium",
		"Pb": "lead",
		"Bi": "bismuth",
		"Po": "polonium",
		"At": "astatine",
		"Rn": "radon",
		"Fr": "francium",
		"Ra": "radium",
		"Rf": "rutherfordium",
		"Db": "dubnium",
		"Sg": "seaborgium",
		"Bh": "bohrium",
		"Hs": "hassium",
		"Mt": "meitnerium",
		"Ds": "darmstadtium",
		"Rg": "roentgenium",
		"Cn": "copemicium",
		"La": "lanthanum",
		"Ce": "cerium",
		"Pr": "praseodymium",
		"Nd": "neodymium",
		"Pm": "promethium",
		"Sm": "samarium",
		"Eu": "europium",
		"Gd": "gadolinium",
		"Tb": "terbium",
		"Dy": "dysprosium",
		"Ho": "holmium",
		"Er": "erbium",
		"Tm": "thulium",
		"Yb": "ytterbium",
		"Lu": "lutetium",
		"Ac": "actinium",
		"Th": "thorium",
		"Pa": "protactinium",
		"U":  "uranium",
		"Np": "neptunium",
		"Pu": "plutonium",
		"Am": "americium",
		"Cm": "curium",
		"Bk": "berkelium",
		"Cf": "californium",
		"Es": "einsteinium",
		"Fm": "fermium",
		"Md": "mendelevium",
		"No": "nobelium",
		"Lr": "lawrencium",
		"D":  "deuterium"}
	// basic molecule regex (atoms and numbers only)
	basicMolRe string
	// (AYZ)n molecule like regex
	oneGroupMolRe string
)

type ByLength []string

func (s ByLength) Len() int {
	return len(s)
}
func (s ByLength) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}
func (s ByLength) Less(i, j int) bool {
	return len(s[i]) > len(s[j])
}

func init() {
	sortedAtoms := make([]string, 0, len(atoms))
	for k := range atoms {
		sortedAtoms = append(sortedAtoms, k)
	}
	// the atom must be sorted by decreasing size
	// to match first Cl before C for example
	sort.Sort(ByLength(sortedAtoms))

	// building the basic molecule regex
	// (atom1|atom2|...)([1-9]*)
	var buffer bytes.Buffer
	buffer.WriteString("(")
	for _, a := range sortedAtoms {
		buffer.WriteString(a)
		buffer.WriteString("|")
	}
	// removing the last |
	buffer.Truncate(buffer.Len() - 1)
	buffer.WriteString(")")
	buffer.WriteString("([1-9]+)*")
	basicMolRe = buffer.String()

	// building the one group molecule regex
	buffer.Reset()
	buffer.WriteString("(?:\\(|\\[)")
	buffer.WriteString("((?:[")
	for _, a := range sortedAtoms {
		buffer.WriteString(a)
		buffer.WriteString("|")
	}
	// removing the last |
	buffer.Truncate(buffer.Len() - 1)
	buffer.WriteString("]+[1-9]*)+)")
	buffer.WriteString("(?:\\)|\\])")
	buffer.WriteString("([1-9]*)")
	oneGroupMolRe = buffer.String()
}

// LinearToEmpiricalFormula returns the empirical formula from the linear formula f.
// Method exported to Python.
// example: [(CH3)2SiH]2NH
//          (CH3)2C[C6H2(Br)2OH]2
//export LinearToEmpiricalFormula
func LinearToEmpiricalFormula(txt *C.char) *C.char {

	// Convert C types to Go types
	f := C.GoString(txt)

	var ef string

	s := "-"
	nf := ""

	// Finding the first (XYZ)n match
	reg := regexp.MustCompile(oneGroupMolRe)

	for s != "" {
		s = reg.FindString(f)

		// Counting the atoms and rebuilding the molecule string
		m := oneGroupAtomCount(s)
		ms := "" // molecule string
		for k, v := range m {
			if v == 1 {
				ms = fmt.Sprintf("%s%s", ms, k)
			} else {
				ms = fmt.Sprintf("%s%s%d", ms, k, v)
			}
		}

		// Then replacing the match with the molecule string - nf is for "new f"
		nf = strings.Replace(f, s, ms, 1)
		f = nf
	}

	// Here nf is a sequence of atoms
	// We just need to count them
	for k, v := range basicAtomCount(nf) {
		if v == 1 {
			ef = ef + fmt.Sprintf("%s", k)
		} else {
			ef = ef + fmt.Sprintf("%s%d", k, v)
		}
	}

	return C.CString(ef)
}

// oneGroupAtomCount returns a count of the atoms of the f formula as a map.
// f must be a formula like (XYZ) (XYZ)n or [XYZ] [XYZ]n.
// example:
// (CH3)2 will return "C":2, "H":6
// CH3CH(NO2)CH3 will return "N":1 "O":2
// CH3CH(NO2)(CH3)2 will return "N":1 "O":2 - process only the first match
func oneGroupAtomCount(f string) map[string]int {
	var (
		// the result map
		c   = make(map[string]int)
		r   *regexp.Regexp
		err error
	)
	// Looking for non matching molecules.
	if r, err = regexp.Compile(oneGroupMolRe); err != nil {
		fmt.Println("error compiling regex:" + oneGroupMolRe)
		return nil
	}
	if !r.MatchString(f) {
		return nil
	}

	// sl is a list of 3 elements like
	// [[(CH3Na6CCl5H)2 CH3Na6CCl5H 2]]
	sl := r.FindAllStringSubmatch(f, -1)
	basicMol := sl[0][1]
	multiplier, _ := strconv.Atoi(sl[0][2])

	// if there is no multiplier
	if multiplier == 0 {
		multiplier = 1
	}

	// counting the atoms
	aCount := basicAtomCount(basicMol)
	for at, nb := range aCount {
		c[at] = nb * multiplier
	}

	return c
}

// basicAtomCount returns a count of the atoms of the f formula as a map.
// f must be a basic formula with only atoms and numbers.
// example:
// C6H5COC6H4CO2H will return "C1":4, "H":10, "O":3
// CH3CH(NO2)CH3 will return Nil, parenthesis are not allowed
func basicAtomCount(f string) map[string]int {
	var (
		// the result map
		c   = make(map[string]int)
		r   *regexp.Regexp
		err error
	)
	// Looking for non matching molecules.
	if r, err = regexp.CompilePOSIX(basicMolRe); err != nil {
		fmt.Println("error compiling regex:" + basicMolRe)
		return nil
	}
	if !r.MatchString(f) {
		return nil
	}

	// sl is a slice like [[Na Na ] [Cl Cl ] [C2 C 2] [Cl3 Cl 3]]
	// for f = NaClC2Cl3
	// [ matchingString capture1 capture2 ]
	// capture1 is the atom
	// capture2 is the its number
	sl := r.FindAllStringSubmatch(f, -1)
	for _, i := range sl {
		atom := i[1]
		var nbAtom int
		if i[2] != "" {
			nbAtom, err = strconv.Atoi(i[2])
			if err != nil {
				return nil
			}
		} else {
			nbAtom = 1
		}
		if _, ok := c[atom]; ok {
			c[atom] = c[atom] + nbAtom
		} else {
			c[atom] = nbAtom
		}
	}
	return c
}

func main() {
}
