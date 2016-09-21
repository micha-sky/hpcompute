require 'csv'

namespace :app do
  namespace :read do
    desc 'Read outloc file'
    task :outlocs, [:river, :file_name] => :environment do |t, args|
      file_name = args.file_name
      river = args.river
      points = []

      times = []
      thalwegs = []
      depths = []
      discharges = []
      file = CSV.open(file_name, {:headers => true, :col_sep => "\t", :skip_blanks => true})
      file.drop(3).each do |row|
        times << row[0].split(' ')[0]
        thalwegs << row[0].split(' ')[3]
        depths << row[0].split(' ')[5]
        discharges << row[0].split(' ')[6]
      end
    end
  end
end