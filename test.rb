
require "./chebyshev.rb"

func = Proc.new { |x| 0.5 * x * x  +  1000.0 * x }
deriv = Proc.new { |x| x + 1000.0 }
integ = Proc.new { |x| x * x * x / 6.0 + 500.0 * x * x }

ch = Chebyshev.new(min: 10.0, max: 100.0, &func)
chint = ch.integral
chder = ch.derivative

ch2 = Chebyshev.new(min: 10.0, max: 100.0, &integ)

t = 10.0
puts "t\tf'(t)\tch'(t)\tdelta"

deriv_min = Float::INFINITY
deriv_max = 0.0

integ_min = Float::INFINITY
integ_max = 0.0

while t <= 100.0
  #puts [t, deriv.call(t), chder.at(t), deriv.call(t) - chder.at(t)].join("\t")
  puts [t, integ.call(t), chint.at(t), integ.call(t) - chint.at(t)].join("\t")
  
  d = deriv.call(t) - chder.at(t)
  deriv_min = d.abs if d.abs < deriv_min
  deriv_max = d.abs if d.abs > deriv_max

  d = integ.call(t) - chint.at(t) - integ.call(10.0)
  integ_min = d.abs if d.abs < integ_min && d != 0.0
  integ_max = d.abs if d.abs > integ_max
  
  t += 0.1
end

puts [deriv_min, deriv_max, integ_min, integ_max].join("\t")

binding.pry
