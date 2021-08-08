package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
)

const (
	gapPenalty = -4
	matrixSize = 24
	inf        = 10000000
)

var (
	proteinMap = map[string]int{
		"A": 0,
		"R": 1,
		"N": 2,
		"D": 3,
		"C": 4,
		"Q": 5,
		"E": 6,
		"G": 7,
		"H": 8,
		"I": 9,
		"L": 10,
		"K": 11,
		"M": 12,
		"F": 13,
		"P": 14,
		"S": 15,
		"T": 16,
		"W": 17,
		"Y": 18,
		"V": 19,
		"B": 20,
		"Z": 21,
		"X": 22,
		"-": 23,
	}
)

func main() {
	// read blosum62.txt
	blosum, err := readBlosum("blosum62.txt")
	if err != nil {
		log.Fatal(err)
	}

	var sc = bufio.NewScanner(os.Stdin)
	aminos := make([]string, 0)
	for sc.Scan() {
		s := sc.Text()
		aminos = append(aminos, s)
	}
	answer, err := star(blosum, aminos)
	if err != nil {
		log.Fatal(err)
	}
	for _, v := range answer {
		fmt.Println(v)
	}
}

func readBlosum(filename string) ([][]int, error) {
	blosum := make([][]int, matrixSize)
	f, err := os.Open(filename)
	if err != nil {
		return nil, err
	}

	r := csv.NewReader(f)
	i := 0
	for {
		record, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if len(record) != matrixSize {
			return nil, fmt.Errorf("error not match matrix Size: %d", len(record))
		}
		blosum[i] = make([]int, len(record))
		for j, str := range record {
			score, err := strconv.Atoi(strings.TrimSpace(str))
			if err != nil {
				return nil, err
			}
			blosum[i][j] = score
		}
		i++
	}
	if len(blosum) != matrixSize {
		return nil, fmt.Errorf("error not match matrix Size: %d", len(blosum))
	}
	return blosum, nil
}

func sequenceAlignmentDP(matrix [][]int, xAmino string, yAmino string) (string, string, int, error) {
	dp := make([][]int, len(xAmino)+1)
	for i := 0; i < len(xAmino)+1; i++ {
		dp[i] = make([]int, len(yAmino)+1)
	}
	for i, x := range xAmino {
		for j, y := range yAmino {
			score, err := getScore(matrix, string(x), string(y))
			if err != nil {
				return "", "", 0, err
			}
			right := dp[i+1][j]
			down := dp[i][j+1]
			if i+1 != len(xAmino) {
				right += gapPenalty
			}
			if j+1 != len(yAmino) {
				down += gapPenalty
			}
			rightDown := dp[i][j] + score
			dp[i+1][j+1] = max3(right, down, rightDown)
		}
	}
	// trace back
	xAlign := ""
	yAlign := ""

	xNow := len(xAmino)
	yNow := len(yAmino)
	for {
		if xNow == 0 && yNow == 0 {
			break
		} else if xNow == 0 {
			xAlign = "-" + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			yNow--
			continue
		} else if yNow == 0 {
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = "-" + yAlign
			xNow--
			continue
		}
		score, err := getScore(matrix, string(xAmino[xNow-1]), string(yAmino[yNow-1]))
		if err != nil {
			return "", "", 0, err
		}
		down := dp[xNow-1][yNow]
		rightDown := dp[xNow-1][yNow-1] + score
		if yNow != len(yAmino) {
			down += gapPenalty
		}
		if dp[xNow][yNow] == rightDown {
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			xNow--
			yNow--
		} else if dp[xNow][yNow] == down {
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = "-" + yAlign
			xNow--
		} else {
			xAlign = "-" + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			yNow--
		}
	}

	return xAlign, yAlign, dp[len(xAmino)][len(yAmino)], nil
}

func star(matrix [][]int, aminos []string) ([]string, error) {
	similaritySums := make([]int, len(aminos))
	for i, rAmino := range aminos {
		for j, cAmino := range aminos[:i] {
			_, _, score, err := sequenceAlignmentDP(matrix, rAmino, cAmino)
			if err != nil {
				return nil, err
			}
			// fmt.Printf("r: %s\nc: %s\n\n", rAlign, cAlign)
			similaritySums[i] += score
			similaritySums[j] += score
		}
	}

	// search max
	_, seqCId := maxSlice(similaritySums)

	// calc
	ansAlignments := make([]string, len(aminos))
	gapList := make([][]int, len(aminos))
	for i, amino := range aminos {
		if i == seqCId {
			ansAlignments[i] = aminos[seqCId]
			continue
		}
		xAlign, yAlign, _, err := sequenceAlignmentDP(matrix, aminos[seqCId], amino)
		// xAlignとaminos[seqCId]を比較してgapの位置を特定する
		// fmt.Printf("r: %s\nc: %s\n\n", xAlign, yAlign)
		if err != nil {
			return nil, err
		}
		ansAlignments[i] = yAlign
		gapList[i], err = searchGaps(xAlign, aminos[seqCId])
	}

	// insert gap
	flag := 0
	nows := make([]int, len(aminos))
	for i := 0; ; i++ {
		// search mins
		minGapPos := inf
		flag = 0
		for j := 0; j < len(aminos); j++ {
			if nows[j] >= len(gapList[j]) {
				continue
			}
			flag = 1
			if gapList[j][nows[j]] < minGapPos {
				minGapPos = gapList[j][nows[j]]
			}
		}
		if flag == 0 {
			break
		}

		// insert gap
		for j := 0; j < len(aminos); j++ {
			if nows[j] < len(gapList[j]) && gapList[j][nows[j]] == minGapPos {
				nows[j]++
				continue
			}
			ok := i + minGapPos
			ansAlignments[j] = ansAlignments[j][:ok] + "-" + ansAlignments[j][ok:]
		}
	}

	return ansAlignments, nil
}

func max3(x, y, z int) int {
	ans := x
	if ans < y {
		ans = y
	}
	if ans < z {
		ans = z
	}
	return ans
}

func maxSlice(x []int) (int, int) {
	var id, value int
	for i, v := range x {
		if v > value {
			value = v
			id = i
		}
	}
	return value, id
}

func maxLength(x []string) (int, int) {
	var id, length int
	for i, v := range x {
		if len(v) > length {
			length = len(v)
			id = i
		}
	}
	return length, id
}

func getProtainId(x string) (int, error) {
	id, ok := proteinMap[x]
	if !ok {
		return 0, fmt.Errorf("No Match Protain: %s", string(x))
	}
	return id, nil
}

func getScore(matrix [][]int, x, y string) (int, error) {
	xid, err := getProtainId(x)
	if err != nil {
		return 0, err
	}
	yid, err := getProtainId(y)
	if err != nil {
		return 0, err
	}
	return matrix[xid][yid], nil
}

func searchGaps(gap, original string) ([]int, error) {
	gapPoint := make([]int, len(gap)-len(original))
	now := 0
	nowGap := 0
	for _, gRune := range gap {
		if now == len(original) {
			gapPoint[nowGap] = now
			nowGap++
			continue
		}
		if string(gRune) != string(original[now]) {
			gapPoint[nowGap] = now
			nowGap++
		} else {
			now++
		}
	}
	if nowGap != len(gap)-len(original) {
		return nil, fmt.Errorf("error: cant find gap")
	}
	return gapPoint, nil
}