use bigdecimal::Num;
use bigdecimal::*;
use std::clone::Clone;
use std::fmt;
use std::ops::{Index, IndexMut, Neg};
use std::slice::{Iter, IterMut};

pub trait NumMatr<T>:
    fmt::Debug + Num + Clone + fmt::Display + Neg<Output = T> + PartialOrd
{
}

impl<T> NumMatr<T> for T where
    T: fmt::Debug + Num + Clone + fmt::Display + Neg<Output = T> + PartialOrd
{
}

pub struct Matrix<T: NumMatr<T> + Signed> {
    data: Box<[T]>,
    rows: usize,
    cols: usize,
}

impl<T: NumMatr<T> + Signed> Matrix<T> {
    pub fn new(rows: usize, cols: usize) -> Self {
        let data = vec![T::zero(); rows * cols];
        Matrix {
            data: data.into_boxed_slice(),
            rows: rows,
            cols: cols,
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
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

    pub fn to_reduced_row_echelon(mut self) -> Self {
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

        // swap rows with bubble sort
        loop {
            let mut swapped = false;
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
                    swapped = true;
                }
            }
            if !swapped {
                break;
            }
        }
        self
    }

    pub fn to_canonical_form(mut self) -> Self {
        for row in 0..self.rows {
            let first_at_row = {
                let mut res = self[(&row, &0)].clone();
                for col in 0..self.cols {
                    if self[(&row, &col)] != T::zero() {
                        res = self[(&row, &col)].clone();
                        break;
                    }
                }
                if res.is_zero() {
                    None
                } else {
                    Some(res)
                }
            };
            if let Some(first_at_row) = first_at_row {
                self.divide_row_by(row, first_at_row);
            }
        }
        self
    }

    pub fn are_zeros_at_maindiag(&self) -> bool {
        for i in 0..std::cmp::min(self.rows, self.cols) {
            if self[(&i, &i)].is_zero() {
                return true;
            }
        }
        return false;
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
        if row1 >= self.rows || row2 >= self.rows {
            panic!();
        }
        for col in 0..self.cols {
            let val1 = self[(&row1, &col)].clone();
            self[(&row1, &col)] = self[(&row2, &col)].clone();
            self[(&row2, &col)] = val1;
        }
    }

    pub fn solve_simple_iterations(&self, accuracy: &T) -> (Vec<T>, usize) {
        if self.cols != self.rows + 1 {
            panic!();
        }
        let mut res: Vec<T> = Vec::with_capacity(self.rows);
        for _ in 0..self.rows {
            res.push(T::zero());
        }

        // do iter
        let mut iters = 0;
        while calc_accuracy(&res, self) > accuracy.clone() && iters < 10000 {
            let newres = res.clone();
            for i in 0..self.rows {
                let mut sum = self[(&i, &(self.cols - 1))].clone();
                for col in 0..(self.cols - 1) {
                    if i != col {
                        sum = sum - self[(&i, &col)].clone() * newres[col].clone();
                    }
                }
                res[i] = sum / self[(&i, &i)].clone();
            }
            iters += 1;
            //println!("{:?}",res);
        }

        (res, iters)
    }

    pub fn solve_zeidel_iterations(&self, accuracy: &T) -> (Vec<T>, usize) {
        if self.cols != self.rows + 1 {
            panic!();
        }
        let mut res: Vec<T> = Vec::with_capacity(self.rows);
        for _ in 0..self.rows {
            res.push(T::zero());
        }

        // do iter
        let mut iters = 0;
        while calc_accuracy(&res, self) > accuracy.clone() && iters < 10000 {
            for i in 0..self.rows {
                let mut sum = self[(&i, &(self.cols - 1))].clone();
                for col in 0..(self.cols - 1) {
                    if i != col {
                        sum = sum - self[(&i, &col)].clone() * res[col].clone();
                    }
                }
                res[i] = sum / self[(&i, &i)].clone();
            }
            iters += 1;
            //println!("accuracy: {:.5}",calc_accuracy(&res, self));
        }

        (res, iters)
    }
}

impl<T: NumMatr<T> + Signed> Index<(&usize, &usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, coords: (&usize, &usize)) -> &Self::Output {
        let row = coords.0;
        let col = coords.1;
        if row >= &self.rows || col >= &self.cols {
            panic!(
                "error: row {}, col {} rows {}, cols {}",
                row, col, self.rows, self.cols
            );
        }
        &self.data[row * self.cols + col]
    }
}

impl<T: NumMatr<T> + Signed> IndexMut<(&usize, &usize)> for Matrix<T> {
    fn index_mut(&mut self, coords: (&usize, &usize)) -> &mut Self::Output {
        let row = coords.0;
        let col = coords.1;
        &mut self.data[row * self.cols + col]
    }
}

impl<T: NumMatr<T> + Signed> fmt::Display for Matrix<T> {
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

impl<T: NumMatr<T> + Signed> Clone for Matrix<T> {
    fn clone(&self) -> Self {
        Matrix {
            cols: self.cols,
            rows: self.rows,
            data: self.data.clone(),
        }
    }
}

pub fn calc_accuracy<T: NumMatr<T> + Signed>(res: &Vec<T>, matr: &Matrix<T>) -> T {
    let mut error = T::zero();
    for row in 0..matr.rows() {
        //println!("row {}", row);
        let result = {
            let mut sum = T::zero();
            for col in 0..(matr.cols() - 1) {
                //println!("col {} {}", col,matr_orig.cols()-1);
                sum = sum + matr[(&row, &col)].clone() * res[col].clone();
            }
            sum
        };
        error = error + (result - matr[(&row, &(matr.cols() - 1))].clone()).abs();
    }
    error
}
