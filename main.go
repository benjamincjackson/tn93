package main

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"runtime"
	"strconv"
	"strings"
	"sync"

	"github.com/spf13/cobra"
)

type ExtendedEncodedFastaRecord struct {
	ID          string
	Description string
	Seq         []byte
	Idx         int
	Count_A     int
	Count_T     int
	Count_G     int
	Count_C     int
}

func ReadEncodeAlignmentToList(f io.Reader) ([]ExtendedEncodedFastaRecord, error) {

	var err error

	records := make([]ExtendedEncodedFastaRecord, 0)

	coding := MakeEncodingArray()

	s := bufio.NewScanner(f)

	first := true

	var id string
	var description string
	var seqBuffer []byte
	var line []byte
	var nuc byte
	var width int
	var counting [256]int

	counter := 0

	for s.Scan() {
		line = s.Bytes()

		if first {

			if len(line) == 0 || line[0] != '>' {
				return []ExtendedEncodedFastaRecord{}, errors.New("badly formatted fasta file")
			}

			description = string(line[1:])
			id = strings.Fields(description)[0]

			first = false

		} else if line[0] == '>' {

			if counter == 0 {
				width = len(seqBuffer)
			} else if len(seqBuffer) != width {
				return []ExtendedEncodedFastaRecord{}, errors.New("different length sequences in input file: is this an alignment?")
			}

			fr := ExtendedEncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
			fr.Count_A = counting[136]
			fr.Count_T = counting[24]
			fr.Count_G = counting[72]
			fr.Count_C = counting[40]
			records = append(records, fr)
			counter++

			description = string(line[1:])
			id = strings.Fields(description)[0]
			seqBuffer = make([]byte, 0)
			for i := range counting {
				counting[i] = 0
			}

		} else {
			encodedLine := make([]byte, len(line))
			for i := range line {
				nuc = coding[line[i]]
				if nuc == 0 {
					return []ExtendedEncodedFastaRecord{}, fmt.Errorf("invalid nucleotide in fasta file (%s)", string(line[i]))
				}
				encodedLine[i] = nuc
				counting[nuc]++
			}
			seqBuffer = append(seqBuffer, encodedLine...)
		}
	}

	if len(seqBuffer) > 0 {
		if counter > 0 && len(seqBuffer) != width {
			return []ExtendedEncodedFastaRecord{}, errors.New("different length sequences in input file: is this an alignment?")
		}
		fr := ExtendedEncodedFastaRecord{ID: id, Description: description, Seq: seqBuffer, Idx: counter}
		fr.Count_A = counting[136]
		fr.Count_T = counting[24]
		fr.Count_G = counting[72]
		fr.Count_C = counting[40]
		records = append(records, fr)
		counter++
	}

	if counter == 0 {
		return []ExtendedEncodedFastaRecord{}, errors.New("empty fasta file")
	}

	err = s.Err()
	if err != nil {
		return []ExtendedEncodedFastaRecord{}, err
	}

	return records, nil
}

func MakeEncodingArray() [256]byte {
	var byteArray [256]byte

	for i := 0; i < 256; i++ {
		byteArray[i] = 0
	}

	byteArray['A'] = 136
	byteArray['a'] = 136
	byteArray['G'] = 72
	byteArray['g'] = 72
	byteArray['C'] = 40
	byteArray['c'] = 40
	byteArray['T'] = 24
	byteArray['t'] = 24
	byteArray['R'] = 192
	byteArray['r'] = 192
	byteArray['M'] = 160
	byteArray['m'] = 160
	byteArray['W'] = 144
	byteArray['w'] = 144
	byteArray['S'] = 96
	byteArray['s'] = 96
	byteArray['K'] = 80
	byteArray['k'] = 80
	byteArray['Y'] = 48
	byteArray['y'] = 48
	byteArray['V'] = 224
	byteArray['v'] = 224
	byteArray['H'] = 176
	byteArray['h'] = 176
	byteArray['D'] = 208
	byteArray['d'] = 208
	byteArray['B'] = 112
	byteArray['b'] = 112
	byteArray['N'] = 240
	byteArray['n'] = 240
	byteArray['-'] = 244
	byteArray['?'] = 242

	return byteArray
}

// See equation (7) in Tamura K, Nei M. Estimation of the number of nucleotide substitutions in the control region of mitochondrial DNA in humans and chimpanzees.
// Mol Biol Evol. 1993 May;10(3):512-26. doi: 10.1093/oxfordjournals.molbev.a040023. PMID: 8336541.
// See also https://github.com/cran/ape/blob/c2fd899f66d6493a80484033772a3418e5d706a4/src/dist_dna.c
func tn93Distance(query, target ExtendedEncodedFastaRecord) float64 {

	// Total ATGC length of the two sequences
	L := float64(target.Count_A + target.Count_C + target.Count_G + target.Count_T + query.Count_A + query.Count_C + query.Count_G + query.Count_T)

	// estimates of the equilibrium base contents from the pair's sequence data
	g_A := float64(target.Count_A+query.Count_A) / L
	g_C := float64(target.Count_C+query.Count_C) / L
	g_G := float64(target.Count_G+query.Count_G) / L
	g_T := float64(target.Count_T+query.Count_T) / L

	g_R := float64(target.Count_A+query.Count_A+target.Count_G+query.Count_G) / L
	g_Y := float64(target.Count_C+query.Count_C+target.Count_T+query.Count_T) / L

	// tidies up the equations a bit, after ape
	k1 := 2.0 * g_A * g_G / g_R
	k2 := 2.0 * g_T * g_C / g_Y
	k3 := 2.0 * (g_R*g_Y - g_A*g_G*g_Y/g_R - g_T*g_C*g_R/g_Y)

	count_P1 := 0 // count of transitional differences between purines (A ⇄ G)
	count_P2 := 0 // count of transitional differences between pyramidines (C ⇄ T)

	count_d := 0 // total number of differences
	count_L := 0 // total length of resolved comparison

	// calculate the three types of change from the pairwise comparison
	for i, tNuc := range target.Seq {
		if query.Seq[i]&8 == 8 && query.Seq[i] == tNuc {
			count_L++
		} else if (query.Seq[i]&tNuc) < 16 && query.Seq[i]&8 == 8 && tNuc&8 == 8 { // are the bases different (and known certainly)
			count_d++
			count_L++
			if (query.Seq[i] | tNuc) == 200 { // 1 if one of the bases is adenine and the other one is guanine, 0 otherwise
				count_P1++
			} else if (query.Seq[i] | tNuc) == 56 { // 1 if one of the bases is cytosine and the other one is thymine, 0 otherwise
				count_P2++
			}
		}
	}

	// estimated rates from this pairwise comparison
	P1 := float64(count_P1) / float64(count_L)                   // rate of changes which are transitional differences between purines (A ⇄ G)
	P2 := float64(count_P2) / float64(count_L)                   // rate of changes which are transitional differences between pyramidines (C ⇄ T)
	Q := float64(count_d-(count_P1+count_P2)) / float64(count_L) // rate of changes which are transversional differences  (A ⇄ C || A ⇄ T || C ⇄ A || C ⇄ G) (i.e. everything else)

	// tidies up the equations a bit, after ape
	w1 := 1.0 - P1/k1 - Q/(2*g_R)
	w2 := 1.0 - P2/k2 - Q/(2*g_Y)
	w3 := 1.0 - Q/(2*g_R*g_Y)

	// tn93 distance:
	d := -k1*math.Log(w1) - k2*math.Log(w2) - k3*math.Log(w3)

	if d == 0.0 {
		d = 0.0
	}

	return d
}

type Pair struct {
	query  ExtendedEncodedFastaRecord
	target ExtendedEncodedFastaRecord
	idx    int
}

type Result struct {
	queryID  string
	targetID string
	d        float64
	idx      int
}

func processPairs(cPair chan Pair, cResults chan Result) {
	for pair := range cPair {
		d := tn93Distance(pair.query, pair.target)
		cResults <- Result{queryID: pair.query.ID, targetID: pair.target.ID, d: d, idx: pair.idx}
	}
}

func printOutput(cResult chan Result, cWriteDone chan bool) {

	outputMap := make(map[int]Result)
	counter := 0

	os.Stdout.WriteString("sequence1\tsequence2\tdistance\n")

	var err error

	for result := range cResult {
		outputMap[result.idx] = result

		if res, ok := outputMap[counter]; ok {
			_, err = os.Stdout.WriteString(res.queryID + "\t" + res.targetID + "\t" + strconv.FormatFloat(res.d, 'f', 10, 64) + "\n")
			if err != nil {
				panic(err)
			}
			delete(outputMap, counter)
			counter++
		} else {
			continue
		}
	}

	for n := 1; n > 0; {
		if len(outputMap) == 0 {
			n--
			break
		}
		res := outputMap[counter]
		_, err = os.Stdout.WriteString(res.queryID + "\t" + res.targetID + "\t" + strconv.FormatFloat(res.d, 'f', 10, 64) + "\n")
		if err != nil {
			panic(err)
		}
		delete(outputMap, counter)
		counter++
	}

	cWriteDone <- true
}

func tn93(infile string, t int) error {

	var threads int
	if t <= 0 {
		threads = runtime.NumCPU()
	} else {
		threads = t
	}

	f, err := os.Open(infile)
	if err != nil {
		return err
	}

	sequences, err := ReadEncodeAlignmentToList(f)
	if err != nil {
		return err
	}

	cPair := make(chan Pair, threads)
	cResults := make(chan Result, threads)
	cWriteDone := make(chan bool)

	go printOutput(cResults, cWriteDone)

	idx_counter := 0
	go func() {
		for i := 0; i < len(sequences)-1; i++ {
			for j := i + 1; j < len(sequences); j++ {
				cPair <- Pair{query: sequences[i], target: sequences[j], idx: idx_counter}
				idx_counter++
			}
		}

		close(cPair)
	}()

	var wgPairs sync.WaitGroup
	wgPairs.Add(threads)

	for n := 0; n < threads; n++ {
		go func() {
			processPairs(cPair, cResults)
			wgPairs.Done()
		}()
	}

	wgPairs.Wait()
	close(cResults)

	_ = <-cWriteDone

	return nil
}

var flagThreads string

func init() {
	mainCmd.Flags().StringVarP(&flagThreads, "threads", "t", "-1", "Number of threads to use (Default: number of available CPUs)")
}

var mainCmd = &cobra.Command{
	Use:   "",
	Short: "",
	Long:  ``,
	RunE: func(cmd *cobra.Command, args []string) error {
		if len(args) != 1 {
			fmt.Println("Usage: ./tn93 -t8 in.fasta > out.tsv")
			os.Exit(1)
		}
		infile := args[0]
		threads, err := strconv.Atoi(flagThreads)
		if err != nil {
			fmt.Println("Usage: ./tn93 -t8 in.fasta > out.tsv")
			return err
		}
		err = tn93(infile, threads)

		return err
	},
}

func usage(*cobra.Command) error {
	fmt.Println("Usage: ./tn93 -t8 in.fasta > out.tsv")
	return nil
}

func main() {
	mainCmd.SetUsageFunc(usage)
	if err := mainCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
