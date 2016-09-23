namespace :app do
  namespace :river do
    desc 'Import river points from kml file'
    task :import, [:river, :file_name] => :environment do |t, args|
      file_name = args.file_name
      river = args.river
      points = []

      doc = Nokogiri::XML(File.open(file_name))

      doc.css('Placemark').each do |placemark|
        # name = placemark.css('Point')
        coordinates = placemark.at_css('coordinates')
        branch = placemark.css('SimpleData')[0].text.to_i
        point = placemark.css('SimpleData')[1].text.to_i
        puts branch

        if  coordinates
          coordinates.text.split(' ').each do |coordinate|
            (lon,lat) = coordinate.split(',')
            point = RiverPoint.new(:river => river, :branch => branch, :point => point,
            :latitude => lat.to_f, :longitude => lon.to_f)
            point.save
            points << point
            print "#{branch},#{point}\n"
          end
        end
      end

    end
  end
end