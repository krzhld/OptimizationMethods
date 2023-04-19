import first_order as first
import second_order as second
import functions as f

if __name__ == '__main__':
    # first.method_of_steepest_descent(f.f, f.grad_f, 0.1)
    second.newton_method(f.f, f.grad_f, f.hess_f, 1e-3)
