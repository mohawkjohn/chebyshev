require 'pry'

class Chebyshev

  # Given a function f, compute the n coefficients of a Chebyshev fit.
  #
  # This code is modified from Numerical Recipes, and I'm not in charge of its license.
  #
  def initialize min: 0.0, max: 1.0, degree: 50, coefficients: Array.new(degree), &f
    @min    = min
    @max    = max
    @degree = degree
    ary     = Array.new(degree)
    @coefficients = coefficients

    if @coefficients.first.nil?
      bma = 0.5 * (max - min)
      bpa = 0.5 * (max + min)

      degree.times do |k|
        y = Math.cos(Math::PI * (k + 0.5) / degree)
        ary[k] = f.call(y * bma + bpa)
      end

      fac = 2.0 / degree
      degree.times do |j|
        sum = 0.0
        degree.times do |k|
          sum += ary[k] * Math.cos(Math::PI * j * (k + 0.5) / degree)
        end

        @coefficients[j] = fac * sum
      end
    end
  end

  attr_reader :min, :max, :degree, :coefficients


  # Evaluate at x
  def at x
    d = 0.0
    dd = 0.0

    raise(RangeError, "Expected between #{min} and #{max} inclusive") if (x - min) * (x - max) > 0.0

    y  = (2.0 * x - min - max) / (max - min).to_f
    y2 = 2.0 * y

    (1..(degree-1)).to_a.reverse.each do |j|
      sv = d
      d = y2 * d - dd + @coefficients[j]
      dd = sv
    end

    y * d - dd + 0.5 * @coefficients[0]
  end


  # Compute the coefficients of the derivative of this Chebyshev approximation.
  def derivative
    c = Array.new(degree)
    c[degree-1] = 0.0
    c[degree-2] = 2.0 * (degree-1) * @coefficients[degree-1]

    (0..(degree-3)).to_a.reverse.each do |j|
      c[j] = c[j+2] + 2 * (j+1) * @coefficients[j+1]
    end

    con = 2.0 / (max - min).to_f
    degree.times do |j|
      c[j] *= con
    end

    Chebyshev.new(coefficients: c, min: min, max: max, degree: degree)
  end


  # Compute the coefficients of the integral of this Chebyshev approximation.
  def integral
    sum = 0.0
    fac = 1.0
    con = 0.25 * (max - min)

    c = Array.new(degree)
    
    (1..(degree-2)).to_a.each do |j|
      c[j]        =  con * (@coefficients[j-1] - @coefficients[j+1]) / j.to_f
      sum        +=  fac * c[j]
      fac         = -fac
    end

    c[degree-1]   = con * @coefficients[degree-2] / (degree-1).to_f
    sum          += fac * c[degree-1]
    c[0]          = 2.0 * sum

    Chebyshev.new(coefficients: c, min: min, max: max, degree: degree)
  end
  
end
