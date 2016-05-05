#!/usr/bin/env ruby
require 'json'

class Coordinate
  attr_accessor :x
  attr_accessor :y
  attr_accessor :z
  def initialize(x, y, z)
    @x = x;
    @y = y;
    @z = z;
  end
end

class Domain
  attr_reader :origin
  attr_reader :x_size
  attr_reader :y_size
  attr_reader :z_size
  attr_reader :volume
  def initialize(x_size, y_size, z_size)
    @origin = Coordinate.new(0, 0, 0)
    @x_size = x_size;
    @y_size = y_size;
    @z_size = z_size;
    @volume = x_size * y_size * z_size;
  end

  def each
    for x in @origin.x...@origin.x+@x_size
      for y in @origin.y...@origin.y+@y_size
        for z in @origin.z...@origin.z+@z_size
          yield x, y, z
        end
      end
    end
  end
end

class Subdomain < Domain
  def initialize(x_size, y_size, z_size, x, y, z)
    super(x_size, y_size, z_size)
    @origin = Coordinate.new(x, y, z);
  end
end

class Separator < Domain
  def initialize(x_size, y_size, z_size, x, y, z)
    super(x_size, y_size, z_size)
    @origin = Coordinate.new(x, y, z)
  end
end

# Divide a domain to two Subdomains and one Separator
def bipart(domain)
  # domain must be a Domain and not a Separator
  if (domain.is_a? Domain) && !(domain.is_a? Separator)
    o = domain.origin
    areas = { x: domain.y_size*domain.z_size, y: domain.x_size*domain.z_size, z: domain.x_size*domain.y_size }
    case areas.min_by(&:last)[0]
    when :x
      lx = (domain.x_size - 1) / 2
      sx = 1
      rx = (domain.x_size - 1) - lx
      ld = Subdomain.new(lx, domain.y_size, domain.z_size, o.x         , o.y, o.z)
      sp = Separator.new( 1, domain.y_size, domain.z_size, o.x + lx    , o.y, o.z)
      rd = Subdomain.new(rx, domain.y_size, domain.z_size, o.x + lx + 1, o.y, o.z)
      return ld, rd, sp
    when :y
      ly = (domain.y_size - 1) / 2
      sy = 1
      ry = (domain.y_size - 1) - ly
      ld = Subdomain.new(domain.x_size, ly, domain.z_size, o.x, o.y         , o.z)
      sp = Separator.new(domain.x_size,  1, domain.z_size, o.x, o.y + ly    , o.z)
      rd = Subdomain.new(domain.x_size, ry, domain.z_size, o.x, o.y + ly + 1, o.z)
      return ld, rd, sp
    when :z
      lz = (domain.z_size - 1) / 2
      sz = 1
      rz = (domain.z_size - 1) - lz
      ld = Subdomain.new(domain.x_size, domain.y_size, lz, o.x, o.y, o.z         )
      sp = Separator.new(domain.x_size, domain.y_size, 1 , o.x, o.y, o.z + lz    )
      rd = Subdomain.new(domain.x_size, domain.y_size, rz, o.x, o.y, o.z + lz + 1)
      return ld ,rd, sp
    end
  end
  return nil, nil, domain
end

# Accept a Domain as an argument and return an Array of subdomains
def recursive_bipart(domain, depth)
  sizes   = []
  domains = [ domain ]
  # Perform bipart until require depth
  lower_bound = 0
  upper_bound = 0
  for l in 0 ... depth
    lower_bound.upto upper_bound do |i|
      left, right, separator = bipart(domains[i])
      sizes[i] = separator.volume
      domains[i] = separator
      domains[2 * i + 1] = left
      domains[2 * i + 2] = right
    end
    lower_bound = 2 * lower_bound + 1
    upper_bound = 2 * upper_bound + 2
  end

  # Update the size of each leaf domains
  length = 2 ** depth
  sizes[lower_bound, length] = domains[lower_bound, length].map(&:volume);

  # Continue to perform bipart for fill reduce reorder
  until domains[lower_bound, length].any? { |d| d.volume < 3 }
    lower_bound.upto upper_bound do |i|
      left, right, separator = bipart(domains[i])
      domains[i] = separator
      domains[2 * i + 1] = left
      domains[2 * i + 2] = right
    end
    # Update bounds
    length = 2 * length
    lower_bound = 2 * lower_bound + 1
    upper_bound = 2 * upper_bound + 2
  end

  # Return a tree of domains & sizes
  return domains, sizes
end

def generate_map(x_size, y_size, z_size, subdomains, depth)
  index = 0
  map = []
  lower_bound = 2 ** depth - 1
  upper_bound = 2 * lower_bound
  factor = (subdomains.length + 1) / (2 ** (depth + 1))
  # Go through each subdomains
  lower_bound.upto upper_bound do |i|
    k = factor
    while k > 0
      for j in (k * i + k - 1) .. (k * i + k * 2 - 2)
        subdomains[j].each do |x, y, z|
          map[index] = x + y * x_size + z * x_size * y_size
          index = index + 1
        end
      end
      k = k / 2
    end
  end

  # Go through each separators
  (depth - 1).downto 0 do |level|
    # Compute new bounds
    lower_bound = (lower_bound - 1) / 2
    upper_bound = (upper_bound - 2) / 2
    lower_bound.upto upper_bound do |i|
      subdomains[i].each do |x, y, z|
        map[index] = x + y * x_size + z * x_size * y_size
        index = index + 1
      end
    end
  end

  return map
end

################################################################################
# Start of this script
################################################################################

# Read parameters in config file
params = JSON.parse(File.read("config.json"))

# Dimensions
x_size = params["x"]
y_size = params["y"]
z_size = params["z"]

depth = params["depth"]

# Output filename
perm_file = params["perm_file"]
part_file = params["part_file"]

if x_size * y_size * z_size > 0
  # Create initial mesh grid
  mesh = Domain.new(x_size, y_size, z_size)

  # Perform recursive bipart (Nested dissection)
  subdomains, sizes = recursive_bipart(mesh, depth)
  
  # Generate reorder mapping
  perm = generate_map(x_size, y_size, z_size, subdomains, depth)

  # Write permutation into perm_file
  File.open(perm_file, 'w') { |file|
    perm.each do |i|
      file.write("#{i}\n")
    end
  }

  # Write partition sizes into perm_file
  File.open(part_file, 'w') { |file|
    sizes.each do |i|
      file.write("#{i}\n")
    end
  }
end
