require File.expand_path('../boot', __FILE__)

require 'rails/all'

# Require the gems listed in Gemfile, including any gems
# you've limited to :test, :development, or :production.
Bundler.require(*Rails.groups)

module Compute
  class Application < Rails::Application
    # Settings in config/environments/* take precedence over those specified here.
    # Application configuration should go into files in config/initializers
    # -- all .rb files in that directory are automatically loaded.

    # Set Time.zone default to the specified zone and make Active Record auto-convert to this zone.
    # Run "rake -D time" for a list of tasks for finding time zone names. Default is UTC.
    # config.time_zone = 'Central Time (US & Canada)'

    # The default locale is :en and all translations from config/locales/*.rb,yml are auto loaded.
    # config.i18n.load_path += Dir[Rails.root.join('my', 'locales', '*.{rb,yml}').to_s]
    # config.i18n.default_locale = :de
    config.sass.preferred_syntax = :sass
    config.autoload_paths += %W(#{config.root}/lib)

    config.assets.paths << "#{Rails.root}/app/assets/videos"
    # Do not swallow errors in after_commit/after_rollback callbacks.
    config.active_record.raise_in_transactional_callbacks = true
  end
end

module Rake
  class Task
    alias_method :origin_invoke, :invoke if method_defined?(:invoke)
    def invoke(*args)
      puts "#{Time.now} STARTED rake task -- #{name} -- #{args.inspect}"
      begin
        result = origin_invoke(*args)
        puts "#{Time.now} ENDED rake task -- #{name}"
        result
      rescue
        puts $!, $@
        puts "#{Time.now} !!! FAILED rake task -- #{name}"
        raise
      end
    end
  end
end
