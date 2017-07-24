package main

import (
	"math/rand"

	"github.com/mingzhi/biogo/seq"
	"github.com/mingzhi/ncbiftp/taxonomy"
)

// SAMPLESIZE define the size of a random sample
const SAMPLESIZE = 100

// CorrResult stores a correlation result.
type CorrResult struct {
	Lag      int
	Mean     float64
	Variance float64
	N        int
	Type     string
}

// CorrResults stores a list of CorrResult with an gene ID.
type CorrResults struct {
	ID      string
	Results []CorrResult
}

// Calculator define a interface for calculating correlations.
type Calculator interface {
	CalcP2(a Alignment, others ...Alignment) (corrResults CorrResults)
}

// CodingCalculator for calculating coding sequences.
type CodingCalculator struct {
	CodingTable *taxonomy.GeneticCode
	MaxCodonLen int
	CodonOffset int
	CodonPos    int
	Synonymous  bool
}

// NewCodingCalculator return a CodingCalculator
func NewCodingCalculator(codingTable *taxonomy.GeneticCode, maxCodonLen, codonOffset, codonPos int, synonymous bool) *CodingCalculator {
	return &CodingCalculator{
		CodingTable: codingTable,
		MaxCodonLen: maxCodonLen,
		CodonOffset: codonOffset,
		CodonPos:    codonPos,
		Synonymous:  synonymous,
	}
}

// CalcP2 calculate P2
func (cc *CodingCalculator) CalcP2(a Alignment, others ...Alignment) CorrResults {
	results := calcP2Coding(a, cc.CodonOffset, cc.CodonPos, cc.MaxCodonLen, cc.CodingTable, cc.Synonymous, calcP2, "P2")
	// p3 := calcP2Coding(a, cc.CodonOffset, cc.CodonPos, cc.MaxCodonLen, cc.CodingTable, cc.Synonymous, calcP3, "P3")
	// p4 := calcP2Coding(a, cc.CodonOffset, cc.CodonPos, cc.MaxCodonLen, cc.CodingTable, cc.Synonymous, calcP4, "P4")
	// results = append(results, p3...)
	// results = append(results, p4...)
	return CorrResults{ID: a.ID, Results: results}
}

func calcP2(codonPairs []CodonPair, codonPos int) (xy float64, n int) {
	nc := doubleCodons(codonPairs, codonPos)
	xy, _, _, n = nc.Cov11()
	return
}

func calcP0(codonPairs []CodonPair, codonPos int) (xy float64, n int) {
	nc := doubleCodons(codonPairs, codonPos)
	xy, _, _, n = nc.Cov00()
	return
}

func calcP3(codonPairs []CodonPair, codonPos int) (xy float64, n int) {
	for n := 0; n < SAMPLESIZE; n++ {
		i := rand.Intn(len(codonPairs))
		j := rand.Intn(len(codonPairs))
		for i == j {
			j = rand.Intn(len(codonPairs))
		}
		k := rand.Intn(len(codonPairs))
		for i == k || j == k {
			k = rand.Intn(len(codonPairs))
		}
		a := codonPairs[i].A[codonPos]
		b := codonPairs[i].B[codonPos]
		c := codonPairs[j].A[codonPos]
		d := codonPairs[k].B[codonPos]
		if a != c && b != d {
			xy++
		}
		n++
	}
	return
}

func calcP4(codonPairs []CodonPair, codonPos int) (xy float64, n int) {
	for n := 0; n < SAMPLESIZE; n++ {
		i := rand.Intn(len(codonPairs))
		j := rand.Intn(len(codonPairs))
		for i == j {
			j = rand.Intn(len(codonPairs))
		}
		k := rand.Intn(len(codonPairs))
		for i == k || j == k {
			k = rand.Intn(len(codonPairs))
		}
		h := rand.Intn(len(codonPairs))
		for i == h || j == h || k == h {
			h = rand.Intn(len(codonPairs))
		}
		a := codonPairs[i].A[codonPos]
		b := codonPairs[h].B[codonPos]
		c := codonPairs[j].A[codonPos]
		d := codonPairs[k].B[codonPos]
		if a != c && b != d {
			xy++
		}
		n++
	}
	return
}

// CalcP is a function type to calculate p2, p3, and p4
type CalcP func(codonPairs []CodonPair, codonPos int) (xy float64, n int)

func calcP2Coding(aln Alignment, codonOffset int, codonPos int, maxCodonLen int, codingTable *taxonomy.GeneticCode, synonymous bool, calcP CalcP, corrType string) (results []CorrResult) {
	codonSequences := [][]Codon{}
	for _, s := range aln.Sequences {
		codons := extractCodons(s, codonOffset)
		codonSequences = append(codonSequences, codons)
	}
	ks := 1.0
	nn := 0
	for l := 0; l < maxCodonLen; l++ {
		totalP2 := 0.0
		totaln := 0
		if l > 0 && ks == 0.0 {
			totalP2 = 0.0
			totaln = nn
		} else {
			for i := 0; i+l < len(codonSequences[0]); i++ {
				codonPairs := []CodonPair{}
				j := i + l
				for _, cc := range codonSequences {
					if i+l < len(cc) {
						codonPairs = append(codonPairs, CodonPair{A: cc[i], B: cc[j]})
					}
				}

				multiCodonPairs := [][]CodonPair{}
				if synonymous {
					multiCodonPairs = synonymousSplit(codonPairs, codingTable)
				} else {
					multiCodonPairs = append(multiCodonPairs, codonPairs)
				}

				for _, codonPairs := range multiCodonPairs {
					if len(codonPairs) >= 2 {
						var codonPositions []int
						if codonPos == -1 {
							codonPositions = []int{0, 1, 2}
						} else {
							codonPositions = append(codonPositions, codonPos)
						}

						for _, codonP := range codonPositions {
							xy, n := calcP(codonPairs, codonP)
							totalP2 += xy
							totaln += n
						}
					}
				}
			}
		}

		if l == 0 {
			ks = totalP2
			nn = totaln
		}

		res1 := CorrResult{
			Lag:  l * 3,
			Mean: totalP2 / float64(totaln),
			N:    totaln,
			Type: corrType,
		}
		results = append(results, res1)
	}

	return
}

func doubleCodons(codonPairs []CodonPair, codonPos int) *NuclCov {
	alphabet := []byte{'A', 'T', 'G', 'C'}
	c := NewNuclCov(alphabet)
	for _, codonPair := range codonPairs {
		a := codonPair.A[codonPos]
		b := codonPair.B[codonPos]
		c.Add(a, b)
	}
	return c
}

// Codon is a byte list of length 3
type Codon []byte

// CodonSequence is a sequence of codons.
type CodonSequence []Codon

// CodonPair is a pair of Codons.
type CodonPair struct {
	A, B Codon
}

// extractCodons return a list of codons from a DNA sequence.
func extractCodons(s seq.Sequence, offset int) (codons []Codon) {
	for i := offset; i+3 <= len(s.Seq); i += 3 {
		c := s.Seq[i:(i + 3)]
		codons = append(codons, c)
	}
	return
}

// synonymousSplit split a list of codon pairs into multiple
// synonymous pairs.
func synonymousSplit(codonPairs []CodonPair, codingTable *taxonomy.GeneticCode) (multiCodonPairs [][]CodonPair) {
	aaList := []string{}
	for _, codonPair := range codonPairs {
		// check gap.
		containsGap := false
		for _, codon := range []Codon{codonPair.A, codonPair.B} {
			for i := 0; i < 3; i++ {
				if codon[i] == '-' || codon[i] == 'N' {
					containsGap = true
					break
				}
			}
		}
		if containsGap {
			continue
		}

		codonA := string(codonPair.A)
		codonB := string(codonPair.B)
		a := codingTable.Table[codonA]
		b := codingTable.Table[codonB]
		ab := string([]byte{a, b})
		index := -1
		for i := 0; i < len(aaList); i++ {
			if aaList[i] == ab {
				index = i
			}
		}
		if index == -1 {
			index = len(aaList)
			aaList = append(aaList, ab)
			multiCodonPairs = append(multiCodonPairs, []CodonPair{})
		}

		multiCodonPairs[index] = append(multiCodonPairs[index], codonPair)
	}

	return
}
