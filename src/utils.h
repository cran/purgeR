#ifndef UTILS_H
#define UTILS_H

template <typename K, typename V>
V map_kv(const std::map <K,V> & m, const K & key, const V & default_value) {
  typename std::map<K,V>::const_iterator it = m.find(key);
  if (it == m.end()) return default_value;
  else return it->second;
}

template<typename T>
int get_index(std::vector<T> v, T i) {
  auto it = find(v.begin(), v.end(), i);
  if (it != v.end()) {
    int index = distance(v.begin(), it);
    return index;
  } else {
    return -1;
  }
}

template<typename T>
inline bool is_in (T i, const std::vector<T>& v) {
  return (std::find(v.begin(), v.end(), i) != v.end());
}

#endif

