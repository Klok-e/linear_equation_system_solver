mod matrix;
use bigdecimal::*;
use matrix::Matrix;
use rand::distributions::{Distribution, Uniform};
use std::env;
use std::fs;

fn main() {
    let mut matr: Matrix<BigDecimal>;

    let path = {
        let mut args: Vec<String> = env::args().collect();
        if args.len() > 1 {
            Some(args.pop().unwrap())
        } else {
            None
        }
    };
    let path = Some("matrix_data.txt");
    if let Some(path) = path {
        let contents = fs::read_to_string(path).expect("Something went wrong reading the file");
        let mut contents = contents
            .split(|x| x == ' ' || x == '\n' || x == ',')
            .map(|x| x.trim())
            .filter(|x| !x.is_empty());

        let columns = contents.next().unwrap().parse::<usize>().unwrap();
        let contents: Vec<f32> = contents.map(|x| x.parse::<f32>().unwrap()).collect();

        println!("{:?}", contents);

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
        matr = Matrix::<BigDecimal>::new(7, 8);

        let mut rng = rand::thread_rng();
        let distr = Uniform::from(-100f32..100f32);

        for elem in matr.iter_mut() {
            *elem = BigDecimal::from(distr.sample(&mut rng));
        }
    }
    println!("{}", matr);
    matr.to_reduced_row_echelon().to_canonical_form();
    println!("{}", matr);
}
