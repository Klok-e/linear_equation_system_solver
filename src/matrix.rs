use bigdecimal::Num;
use std::clone::Clone;
use std::fmt;
use std::ops::{Index, IndexMut, Neg};
use std::slice::{Iter, IterMut};

pub trait NumMatr<T>: Num + Clone + fmt::Display + Neg<Output = T> {}

impl<T> NumMatr<T> for T where T: Num + Clone + fmt::Display + Neg<Output = T> {}

pub struct Matrix<T: NumMatr<T>> {
    data: Box<[T]>,
    rows: usize,
    cols: usize,
}

impl<T: NumMatr<T>> Matrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        let data = vec![T::zero(); rows * cols];
        Matrix {
            data: data.into_boxed_slice(),
            rows: rows,
            cols: cols,
        }
    }

    pub fn with(data: Vec<T>, rows: usize, cols: usize) -> Self {
        if rows * cols != data.len() {
            panic!("wrong data size");
        }
        Matrix {
            data: data.into_boxed_slice(),
            rows: rows,
            cols: cols,
        }
    }

    pub fn iter(&self) -> Iter<T> {
        self.data.iter()
    }

    pub fn iter_mut(&mut self) -> IterMut<T> {
        self.data.iter_mut()
    }

    pub fn to_reduced_row_echelon(&mut self) -> &mut Self {
        // to row echelon form
        for row in 0..self.rows {
            let lead_coef: Option<(T, usize)> = {
                let mut col_ret = None;
                // self.cols - 1 because the last column is of an augmented matrix
                for col in 0..(self.cols - 1) {
                    if self[(&row, &col)] != T::zero() {
                        col_ret = Some((self[(&row, &col)].clone(), col.clone()));
                        break;
                    }
                }
                col_ret
            };
            if let Some((lead_coef, col_ind)) = lead_coef {
                for row_below in (row + 1)..self.rows {
                    let mult = self[(&row_below, &col_ind)].clone();
                    self.add_rows_with_mult_by(row_below, row, lead_coef.clone(), -mult);
                }
            }
        }

        // to reduced row echelon form
        for row in (1..self.rows).rev() {
            let lead_coef: Option<(T, usize)> = {
                let mut col_ret = None;
                // self.cols - 1 because the last column is of an augmented matrix
                for col in 0..(self.cols - 1) {
                    if self[(&row, &col)] != T::zero() {
                        col_ret = Some((self[(&row, &col)].clone(), col.clone()));
                        break;
                    }
                }
                col_ret
            };
            if let Some((lead_coef, col_ind)) = lead_coef {
                for row_above in (0..row).rev() {
                    let mult = self[(&row_above, &col_ind)].clone();
                    self.add_rows_with_mult_by(row_above, row, lead_coef.clone(), -mult);
                }
            }
        }

        // swap rows
        for row in 0..(self.rows - 1) {
            let col_argnonzero = {
                let mut arg = 0;
                for col in 0..self.cols {
                    if self[(&row, &col)] != T::zero() {
                        arg = col;
                        break;
                    }
                }
                arg
            };
            let col_rowbelow_argnonzero = {
                let mut arg = 0;
                for col in 0..self.cols {
                    if self[(&(row + 1), &col)] != T::zero() {
                        arg = col;
                        break;
                    }
                }
                arg
            };
            if col_argnonzero > col_rowbelow_argnonzero {
                self.swap_rows(row, row + 1);
            }
        }
        self
    }

    pub fn to_canonical_form(&mut self) -> &mut Self {
        for row in 0..self.rows {
            let first_at_row = {
                let mut res = T::zero();
                for col in 0..self.cols {
                    if self[(&row, &col)] != T::zero() {
                        res = self[(&row, &col)].clone();
                        break;
                    }
                }
                res
            };
            self.divide_row_by(row, first_at_row);
        }
        self
    }

    fn divide_row_by(&mut self, row: usize, div: T) {
        for col in 0..self.cols {
            self[(&row, &col)] = self[(&row, &col)].clone() / div.clone();
        }
    }

    fn add_rows_with_mult_by(
        &mut self,
        row_to_be_added_to: usize,
        row_to_add: usize,
        row_to_be_added_to_mult_by: T,
        row_to_add_mult_by: T,
    ) {
        for column in 0..self.cols {
            let elem_to_add = self[(&row_to_add, &column)].clone();
            let elem_to_be_added_to = self[(&row_to_be_added_to, &column)].clone();
            self[(&row_to_be_added_to, &column)] = elem_to_be_added_to
                * row_to_be_added_to_mult_by.clone()
                + elem_to_add * row_to_add_mult_by.clone();
        }
    }

    pub fn swap_rows(&mut self, row1: usize, row2: usize) {
        for col in 0..self.cols {
            let val1 = self[(&row1, &col)].clone();
            self[(&row1, &col)] = self[(&row2, &col)].clone();
            self[(&row2, &col)] = val1;
        }
    }
}

impl<T: NumMatr<T>> Index<(&usize, &usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, coords: (&usize, &usize)) -> &Self::Output {
        let row = coords.0;
        let col = coords.1;
        &self.data[row * self.cols + col]
    }
}

impl<T: NumMatr<T>> IndexMut<(&usize, &usize)> for Matrix<T> {
    fn index_mut(&mut self, coords: (&usize, &usize)) -> &mut Self::Output {
        let row = coords.0;
        let col = coords.1;
        &mut self.data[row * self.cols + col]
    }
}

impl<T: NumMatr<T>> fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[")?;
        for (i, elem) in self.iter().enumerate() {
            if i % self.cols == self.cols - 1 {
                write!(
                    f,
                    "{:12.3}]\n{}",
                    elem,
                    if i == self.cols * self.rows - 1 {
                        ""
                    } else {
                        "["
                    }
                )?;
            } else {
                write!(f, "{:12.3}, ", elem)?;
            }
        }
        write!(f, "\n")
    }
}