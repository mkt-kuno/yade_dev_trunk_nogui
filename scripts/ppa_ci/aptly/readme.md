# Specific scripts, applied on yade-dem server

Scripts are placed on the yade-dem server in /home/anton/apt

Steps to add a new distribution:

1. Update [docker-yade](https://gitlab.com/yade-dev/docker-yade) to build the new container.
2. Update [.gitlab-ci](https://gitlab.com/yade-dev/trunk/-/tree/master/.gitlab-ci) to build yade against new image.
3. Update /home/anton/create_repos.sh with a new distribution name.
4. Run the /home/anton/create_repos.sh script.
5. Update the /home/anton/aptly/update_repos_initial.sh script to include the new distribution.
6. Run the /home/anton/aptly/update_repos_initial.sh script.
7. Update the /home/anton/aptly/update_repos_next.sh.