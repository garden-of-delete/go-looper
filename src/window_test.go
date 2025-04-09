package src

import (
	"reflect"
	"testing"
)

func TestFromLinearWindows(t *testing.T) {
	inputSeq := []rune{'G', 'A', 'T'}
	//expected := []Window{}
	test1 := FromLinearWindows(inputSeq, 2)
	test2 := FromLinearWindows(inputSeq, 3)
	test3 := FromLinearWindows(inputSeq, 4)
	expected1 := []Window{{0, 1}, {0, 2}, {1, 2}}
	expected2 := []Window{{0, 2}}
	expected3 := []Window{}

	if !reflect.DeepEqual(test1, expected1) {
		t.Errorf("FromLinearWindows(%v) = %v, want %v", inputSeq, test1, expected1)
	}
	if !reflect.DeepEqual(test2, expected2) {
		t.Errorf("FromLinearWindows(%v) = %v, want %v", inputSeq, test2, expected2)
	}
	if len(test3) != len(expected3) {
		t.Errorf("FromLinearWindows(%v) = %v, want %v", inputSeq, test3, expected3)
	}

}

func TestFromCircularWindows(t *testing.T) {
	inputSeq := []rune{'G', 'A', 'T'}
	inputSeq2 := []rune{'G', 'A', 'T', 'T'}
	test1 := FromCircularWindows(inputSeq, 2)
	test2 := FromCircularWindows(inputSeq2, 2)
	test3 := FromCircularWindows(inputSeq2, 3)
	expected1 := []Window{{1, 0}, {2, 0}, {2, 1}}
	expected2 := []Window{{1, 0}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}}
	expected3 := []Window{{1, 0}, {2, 0}, {2, 1}, {3, 1}, {3, 2}, {2, 3}} // TODO: fix

	if !reflect.DeepEqual(test1, expected1) {
		t.Errorf("FromCircularWindows(%v) = %v, want %v", inputSeq, test1, expected1)
	}
	if !reflect.DeepEqual(test2, expected2) {
		t.Errorf("FromCircularWindows(%v) = %v, want %v", inputSeq2, test2, expected2)
	}
	if !reflect.DeepEqual(test3, expected3) {
		t.Errorf("FromCircularWindows(%v) = %v, want %v", inputSeq2, test3, expected3)
	}
}
