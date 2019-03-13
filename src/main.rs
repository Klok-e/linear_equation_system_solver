mod matrix;
use bigdecimal::*;
use matrix::Matrix;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::SeedableRng;
use rand::rngs::SmallRng;
use std::env;
use std::fs;

/*
First number in a source file is number of columns, rows are inferred
*/

fn main() {
    let path = {
        let mut args: Vec<String> = env::args().collect();
        if args.len() > 1 {
            Some(args.pop().unwrap())
        } else {
            None
        }
    };
    //let path = Some("matrix_data2.txt");

    let matr = {
        let mut matr: Matrix<BigDecimal>;

        if let Some(path) = path {
            let contents = fs::read_to_string(path).expect("Something went wrong reading the file");
            let mut contents = contents
                .split(|x| x == ' ' || x == '\n' || x == ',')
                .map(|x| x.trim())
                .filter(|x| !x.is_empty());

            let columns = contents.next().unwrap().parse::<usize>().unwrap();
            let contents: Vec<f32> = contents.map(|x| x.parse().unwrap()).collect();

            if contents.len() % columns != 0 {
                panic!();
            }
            matr = Matrix::with(
                {
                    let mut vec = Vec::<BigDecimal>::with_capacity(contents.len());
                    for item in contents.iter() {
                        vec.push(BigDecimal::from(*item));
                    }
                    vec
                },
                contents.len() / columns,
                columns,
            );
        } else {
            matr = Matrix::<BigDecimal>::new(6, 7);
            let seed = {
                let mut arr = [0u8; 16];
                arr[0] = 240;
                arr[1] = 140;
                arr[2] = 62;
                arr
            };
            let mut rng = SmallRng::from_seed(seed);
            let distr = Uniform::from(-10i32..10i32);

            for elem in matr.iter_mut() {
                *elem = BigDecimal::from(distr.sample(&mut rng));
            }
        }
        matr
    };
    println!("{}", matr);
    let matr = matr.to_reduced_row_echelon();
    println!("{}", matr);
    let matr = matr.to_canonical_form();
    println!("{}", matr);
}
