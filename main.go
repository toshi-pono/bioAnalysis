package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"strconv"
	"strings"
	"time"

	echo "github.com/labstack/echo/v4"
)

const (
	gapPenalty = -4
	matrixSize = 24
	inf        = 10000000
)

var (
	blosum     [][]int
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

type Tree struct {
	Left  int     `json:"left"`
	Right int     `json:"right"`
	Score float64 `json:"score"`
}

type AminoSequences struct {
	Aminos []string `json:"aminos"`
}

type TreeResponse struct {
	TreeData    []Tree `json:"treedata"`
	ElapsedTime int64  `json:"elapsedtime"`
}

func main() {
	// load blosum62.txt
	var err error
	blosum, err = readBlosum("blosum62.txt")
	if err != nil {
		log.Fatal(err)
	}

	e := echo.New()
	e.GET("/", getIndexHandler)
	e.POST("/tree", postTreeHandler)
	e.Static("/static", "public")
	e.Logger.Fatal(e.Start(":8080"))
}

// getIndexHandler return index.html
func getIndexHandler(c echo.Context) error {
	return c.File("public/index.html")
}

// postTreeHandler create tree
func postTreeHandler(c echo.Context) error {
	param := new(AminoSequences)
	if err := c.Bind(param); err != nil {
		return err
	}

	now := time.Now()
	multiAlignment, err := star(blosum, param.Aminos)
	if err != nil {
		log.Fatal(err)
	}
	tree, err := createTree(multiAlignment)
	if err != nil {
		log.Fatal(err)
	}

	return c.JSON(http.StatusOK, TreeResponse{
		TreeData:    tree,
		ElapsedTime: time.Since(now).Milliseconds(),
	})
}

// readBlosum ファイルからアミノ酸残基間スコアを読み込み
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

// sequenceAlignmentDP 動的計画法によってペア間アライメントとスコアを計算する
func sequenceAlignmentDP(matrix [][]int, xAmino string, yAmino string) (string, string, int, error) {
	dp := make([][]int, len(xAmino)+1)
	for i := 0; i < len(xAmino)+1; i++ {
		dp[i] = make([]int, len(yAmino)+1)
	}
	// DP
	for i, x := range xAmino {
		for j, y := range yAmino {
			score, err := getScore(matrix, string(x), string(y))
			if err != nil {
				return "", "", 0, err
			}
			// 3方向の移動について比較して最大のものを選択
			right := dp[i+1][j]
			down := dp[i][j+1]
			rightDown := dp[i][j] + score

			// Out Gap のペナルティはゼロ
			if i+1 != len(xAmino) {
				right += gapPenalty
			}
			if j+1 != len(yAmino) {
				down += gapPenalty
			}
			dp[i+1][j+1] = max3(right, down, rightDown)
		}
	}
	// trace back
	xAlign := ""
	yAlign := ""

	xNow := len(xAmino)
	yNow := len(yAmino)
	for {
		// 左辺か上辺にたどり着いたとき
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
		// どの方向から来たか大きさを比較して戻る
		score, err := getScore(matrix, string(xAmino[xNow-1]), string(yAmino[yNow-1]))
		if err != nil {
			return "", "", 0, err
		}
		down := dp[xNow-1][yNow]
		rightDown := dp[xNow-1][yNow-1] + score
		// Out Gapペナルティはゼロ
		if yNow != len(yAmino) {
			down += gapPenalty
		}
		if dp[xNow][yNow] == rightDown {
			// 右下
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			xNow--
			yNow--
		} else if dp[xNow][yNow] == down {
			// 下
			xAlign = string(xAmino[xNow-1]) + xAlign
			yAlign = "-" + yAlign
			xNow--
		} else {
			// 右
			xAlign = "-" + xAlign
			yAlign = string(yAmino[yNow-1]) + yAlign
			yNow--
		}
	}

	return xAlign, yAlign, dp[len(xAmino)][len(yAmino)], nil
}

// star スター法によって多重アライメントを計算
func star(matrix [][]int, aminos []string) ([]string, error) {
	// ペア間類似性スコアの計算
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

	// 中心となる配列 seqC の決定
	_, seqCId := maxSlice(similaritySums)

	// 配列 seqC を中心にペア間アライメントを準備
	ansAlignments := make([]string, len(aminos))
	gapList := make([][]int, len(aminos))
	for i, amino := range aminos {
		if i == seqCId {
			ansAlignments[i] = aminos[seqCId]
			continue
		}
		// DP
		xAlign, yAlign, _, err := sequenceAlignmentDP(matrix, aminos[seqCId], amino)
		// fmt.Printf("r: %s\nc: %s\n\n", xAlign, yAlign)
		if err != nil {
			return nil, err
		}

		ansAlignments[i] = yAlign
		// xAlignとaminos[seqCId]を比較してgapの位置を特定する
		gapList[i], err = searchGaps(xAlign, aminos[seqCId])
	}

	// Gapを挿入してアライメントを融合
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

// max3 Return the largest of the three
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

// maxSlice Return the maximum in the slice
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

// maxLength Return the maximum length of the slice
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

// getProtainId Get the number corresponding to the amino acid
func getProtainId(x string) (int, error) {
	id, ok := proteinMap[x]
	if !ok {
		return 0, fmt.Errorf("No Match Protain: %s", string(x))
	}
	return id, nil
}

// getScore Get the similarity score
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

// searchGaps Locate the gap insertion point
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

// pDistance (距離 = 配列の異なる部分の数 / 配列長)
func pDistance(x, y string) (float64, error) {
	if len(x) != len(y) {
		return 0.0, fmt.Errorf("string length is diffrence")
	}
	count := 0
	for i := 0; i < len(x); i++ {
		if x[i] != y[i] {
			count++
		}
	}
	return float64(count) / float64(len(x)), nil
}

// createTree UPGMA法によって進化系統樹を作成する
func createTree(multiAlignment []string) ([]Tree, error) {
	distances := make([][]float64, len(multiAlignment))
	for i := 0; i < len(multiAlignment); i++ {
		distances[i] = make([]float64, len(multiAlignment))
	}

	// 距離行列を作成
	for i := 0; i < len(multiAlignment); i++ {
		for j := 0; j < i; j++ {
			distance, err := pDistance(multiAlignment[i], multiAlignment[j])
			if err != nil {
				return nil, err
			}
			distances[i][j] = distance
			distances[j][i] = distance
		}
	}

	// クラスタに含まれる値の数を保持
	duplications := make([]int, len(multiAlignment))
	for i := 0; i < len(duplications); i++ {
		duplications[i] = 1
	}

	// 木に葉を追加
	nodeIds := make([]int, len(multiAlignment))
	tree := make([]Tree, 0)
	for i := 0; i < len(multiAlignment); i++ {
		tree = append(tree, Tree{
			Left:  -1,
			Right: -1,
			Score: 0,
		})
		nodeIds[i] = i
	}

	// 進化系統樹の作成
	for k := 0; k < len(multiAlignment)-1; k++ {
		// 距離行列の中で，最小のペアを見つける
		minDistance := 1.0
		var minI, minJ int
		for i := 0; i < len(multiAlignment); i++ {
			if duplications[i] == 0 {
				continue
			}
			for j := 0; j < i; j++ {
				if duplications[j] == 0 {
					continue
				}
				if minDistance > distances[i][j] {
					minI = i
					minJ = j
					minDistance = distances[i][j]
				}
			}
		}

		// minJにminIを統合する（距離行列の更新）
		for i := 0; i < len(multiAlignment); i++ {
			if duplications[i] == 0 || i == minI || i == minJ {
				continue
			}
			distance := (distances[minJ][i]*float64(duplications[minJ]) + distances[minI][i]*float64(duplications[minI])) / float64(duplications[minI]+duplications[minJ])
			distances[minJ][i] = distance
			distances[i][minJ] = distance
		}

		// 木に節を追加する
		tree = append(tree, Tree{
			Left:  nodeIds[minJ],
			Right: nodeIds[minI],
			Score: minDistance / 2.0,
		})
		duplications[minJ] += duplications[minI]
		duplications[minI] = 0
		nodeIds[minJ] = len(tree) - 1
	}

	return tree, nil
}
