require 'csv'

namespace :app do
    desc 'Cleans tmp directory'
    task :cleanup, [] => :environment do |t, args|
        Dir.glob('tmp') do |tmp_file|
            if File.mtime(tmp_file) < (Time.now - 7.days)
                File.delete(tmp_file)
                puts "Deleted increment file #{tmp_file}"
            end
        end

    end
end