package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

const (
	gapPenalty = -4
	matrixSize = 24
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
	// var sc = bufio.NewScanner(os.Stdin)
	// read blosum62.txt
	blosum, err := readBlosum("blosum62.txt")
	if err != nil {
		log.Fatal(err)
	}
	xAlign, yAlign, err := sequenceAlignmentDP(blosum, "ARBND", "ARN")
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(xAlign)
	fmt.Println(yAlign)

	// for sc.Scan() {
	// 	s := sc.Text()
	// 	fmt.Println(s)
	// }
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

func sequenceAlignmentDP(matrix [][]int, xAmino string, yAmino string) (string, string, error) {
	dp := make([][]int, len(xAmino)+1)
	for i := 0; i < len(xAmino)+1; i++ {
		dp[i] = make([]int, len(yAmino)+1)
	}
	for i, x := range xAmino {
		for j, y := range yAmino {
			score, err := getScore(matrix, string(x), string(y))
			if err != nil {
				return "", "", err
			}
			dp[i+1][j+1] = max3(dp[i+1][j]+gapPenalty, dp[i][j+1]+gapPenalty, dp[i][j]+score)
		}
	}

	// get back
	size := int(math.Max(float64(len(xAmino)), float64(len(yAmino))))
	xAlign := ""
	yAlign := ""

	xNow := len(xAmino)
	yNow := len(yAmino)
	for i := 0; i < size; i++ {
		score, err := getScore(matrix, string(xAmino[xNow-1]), string(yAmino[yNow-1]))
		if err != nil {
			return "", "", err
		}
		if dp[xNow][yNow] == dp[xNow-1][yNow-1]+score {
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			xNow--
			yNow--
		} else if dp[xNow][yNow] == dp[xNow-1][yNow]+gapPenalty {
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = "-" + yAlign
			xNow--
		} else {
			xAlign = "-" + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			yNow--
		}
	}

	return xAlign, yAlign, nil
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
